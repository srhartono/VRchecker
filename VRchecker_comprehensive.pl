#!/usr/bin/perl

use strict; use warnings; use Getopt::Std; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_i $opt_b $opt_s $opt_w $opt_l $opt_o $opt_H $opt_q $opt_L $opt_t $opt_a);
getopts("vi:b:s:w:l:o:HqL:t:a");

BEGIN {
   my $libPath = dirname(dirname abs_path $0) . '/';
   push(@INC, $libPath);
}

use FAlite;

our $N         ="\e[0m";
our $B         ="\e[1m";
our $BL     ="\e[0;30m";
our $BU     ="\e[0;34m";
our $GN     ="\e[0;32m";
our $YW     ="\e[1;33m";
our $CY     ="\e[0;36m";
our $RD     ="\e[0;31m";
our $PR     ="\e[0;35m";
our $BR     ="\e[0;33m";
our $GR     ="\e[0;37m";
our $WT     ="\e[0m";
our $DGR    ="\e[1;30m";
our $LBU    ="\e[1;34m";
our $LGN    ="\e[1;32m";
our $LCY    ="\e[1;36m";
our $LRD    ="\e[1;31m";
our $LPR    ="\e[1;35m";

my ($fqFile, $VRFile) = ($opt_i, $opt_b);
my $usage = "
Usage: $YW$0$N -i $LCY<fqFile.fastq.gz>$N -b $LPR<VRFile.fa>$N

-t: percent match threshold [75]
-o: output dir [currrent directory]
-q: quiet, don't print anything

Seeding match parameters:
-s: step size [20]
-w: window subtract size [25]
-l: min continuous match length threshold [20]
-L: min continuous match length threshold as of percent VR length [10]
-a: continue searching all VRs for better matches (slower)

- At least 20bp of VR sequence match:
$YW$0$N -i ${LCY}example.fastq.gz$N -b ${LGN}example_VR.fa$N -l ${LPR}20$N

Once seed is found in the read, then quantify the VR's matching score using ${YW}clustalo-1.2.4-Ubuntu-x86_64$N
unless -a, then we only take first VR whose match is above threshold -t [75]

More options & explanation do -H

";

my $usage_long = "

Options for quick dirty sliding-window based non-exact matching (blast takes too long)
-s: step size [20]
-w: window subtract size [25]
-l: min continuous match length threshold [20]
-L: min continuous match length threshold as of percent VR length [10]
-a: continue searching all VRs for better matches (slower)


If the whole VR length doesn't match (200bp) it'll try to reduce the window size by [-w] 
and do sliding window matching with step size of [-s]
It'll stop when finding a match or when window size is less than min length threshold [-l, default 20]

Example: -s step size 20, -w window 25, min match length threshold$LGN 100$N

Initial VR sequence size is 200, try match to the read sequence, if not match then:
window = 200 - 25 = 175 bp -> substring position 0-175 and 20-195 of VR sequence, then try match to read sequence, if not:
window = 175 - 25 = 150 bp -> substring position 0-150, 20-170, and 40-190 of VR sequence
window = 150 - 25 = 125 bp -> substring position 0-125, 20-145, 40 - 165, 60 - 185 of VR sequence
window = 125 - 25 = 100 bp -> substring position 0-100, 20-120, 40-140, 60-160, 80-180, 100-200
window = 100 - 25 =  75 bp -> search stop here coz window length (75bp) is less than min match length threshold -l 100bp

If min match length threshold is same as VR size (e.g. 200bp), then step size/window doesn't matter

Once seed is found in the read, then quantify the VR's matching score using ${YW}clustalo-1.2.4-Ubuntu-x86_64$N
unless -a, then we only take first VR whose match is above threshold -t [75]


";

die $usage unless defined $opt_i and defined $opt_b;
die $usage . "\n" . $usage_long if defined $opt_H;
die "fastq file -i $LCY$fqFile$N doesn't exist!\n" if not -e $fqFile;
die "query file -b $LCY$VRFile$N doesn't exist!\n" if not -e $VRFile;
my $step_size             = defined $opt_s ? $opt_s : 20;
my $window_subtract_size  = defined $opt_w ? $opt_w : 25;
my $minlen_threshold_orig = defined $opt_l ? $opt_l : 20;
my $minlen_threshold_perc = defined $opt_L ? $opt_L : 10;
my $outDir                = defined $opt_o ? $opt_o : "./";
my $perc_match_threshold  = defined $opt_t ? $opt_t : 75;
my $quiet = 1 if defined $opt_q;
#
#" unless @ARGV;
my %VR;
system("mkdir -p $outDir") if not -d $outDir;
$outDir =~ s/\/+$//;
my ($folder1, $fqFilename) = getFilename($fqFile, "folderfull");
my ($folder2, $VRFilename) = getFilename($VRFile, "folderfull");
my $logFile = "$outDir/$fqFilename\_$VRFilename\_l$minlen_threshold_orig\_w$window_subtract_size\_s$step_size.comprehensive.log";
my $outFile = "$outDir/$fqFilename\_$VRFilename\_l$minlen_threshold_orig\_w$window_subtract_size\_s$step_size.comprehensive.tsv";
my $outFileshort = "$outDir/$fqFilename\_$VRFilename\_l$minlen_threshold_orig\_w$window_subtract_size\_s$step_size.comprehensive.short.tsv";
my $VRstatsFile = "$outDir/$fqFilename\_$VRFilename\_l$minlen_threshold_orig\_w$window_subtract_size\_s$step_size.comprehensive.stats";
open (my $outLog, ">", $logFile) or die "Can't write to $logFile: $!\n";
my %out;
LOG($outLog, "$LGN
Parameters:$N
-i fastq file            : $LCY$fqFile$N
-b VR/query fasta file   : $LCY$VRFile$N
-o output dir            : $LCY$outDir$N
-s step size             : $LGN$step_size$N
-w window subtract size  : $LGN$window_subtract_size$N
-l min length threshold  : $LGN$minlen_threshold_orig$N
-t perc match threshold  : $LGN$perc_match_threshold$N

",$opt_q);
LOG($outLog, "\n#VR sequenecs:\n",$opt_q);
open (my $in0, "<", "$VRFile") or die "Cannot read from $VRFile: $!\n\n";
my $fasta = new FAlite($in0);
while (my $entry = $fasta->nextEntry) {
	my ($VRseq) = $entry->seq;
	my ($VRdef) = $entry->def; $VRdef =~ s/^>//;
	$VR{$VRseq} = $VRdef;
	my $VRlen = length($VRseq);
	my $VRshort = $VRseq;
	($VRshort) = $VRseq =~ /^(.{20})/ if length($VRseq) >= 20;
   my $minlen_threshold = $minlen_threshold_orig eq 0 ? length($VRseq) : $minlen_threshold_orig;
   if (defined $opt_L) {
      $minlen_threshold = int($minlen_threshold_perc / 100 * length($VRseq)+0.5);
   }
   print "$LCY$VRdef$N\t$YW$VRshort$N... ($LGN$VRlen$N bp, minlen=$LGN$minlen_threshold$N bp)\n" if not defined $opt_q;
   LOG($outLog, "$VRdef\t$VRseq\n",1);
}
close $in0;
my $outUnk;
LOG($outLog, "\n# matching query sequence to fastq sequence\n",$opt_q);
my $linecount = 0;
my %data;
my $inputFastQ = $fqFile;
open (my $out0, ">", $outFile) or die "Failed to write to $outFile: $!\n";
open (my $out0b, ">", $outFileshort) or die "Failed to write to $outFileshort: $!\n";
my $inFastQ;
if ($inputFastQ =~ /.gz$/) {
	open ($inFastQ, "zcat $inputFastQ|") or die "Cannot read from $inputFastQ: $!\n";
}
else {
	open ($inFastQ, "<", "$inputFastQ") or die "Cannot read from $inputFastQ: $!\n";
}
print $out0 "#read_name\tmatch_start\tmatch_end\tmatch_name\tmatch_number\tmatch_strand\tmatch_seq\tread_seq\tread_qual\n";
while (my $line = <$inFastQ>) {
	chomp($line);
	LOG($outLog, "Done $LGN$linecount$N Reads\n") if $linecount % 100 == 0 and not defined $opt_q;
	$linecount ++;

	my $read_name = $line;
	my ($read_namecheck) = $read_name =~ /^.(.+)$/;

	$line = <$inFastQ>; chomp($line);
	my $read_seq = $line;

	$line = <$inFastQ>; chomp($line);
	my ($pluscheck) = $line =~ /^.(.+)$/;
	if ($line ne "+" and $read_namecheck ne $pluscheck) {
		die "\nfastq corrupted? Line $LGN$linecount$N isn't a + or same name as def ($LCY$read_name$N)\n";
	}
	$line = <$inFastQ>; chomp($line);
	my $read_qual = $line;
	my $VRfinal = "UNKNOWN";
	my $VRseqfinal = "";
	my $VRstrand = ".";
	my ($VRseqbeg, $VRseqend);
	my ($beg, $end) = (-1,-1);
	$data{total} ++;
	my $totalbc = 0;
	my @bc;
	my $score = 0;
	my @bcseq;
	my $bestlength = 0;
	#print "foreach\n" if ($read_name eq "\@6acb2984-5480-415b-96db-da1208fa1a0d");
	my $dontuse;
	foreach my $VRorig (sort keys %VR) {
		my $VRname = $VR{$VRorig};
		next if defined $dontuse->{$VRname};
		my $minlen = length($VRorig);
		my $init = 0;
		my $found = 0;
		my $VRrevorig = revcomp($VRorig);
	   my $minlen_threshold = $minlen_threshold_orig eq 0 ? length($VRorig) : $minlen_threshold_orig;
	   if (defined $opt_L) {
	      $minlen_threshold = int($minlen_threshold_perc / 100 * length($VRorig)+0.5);
	   }
		my ($bcarr, $bcseqarr, $VRstrand2);
		#print "   while $LGN$VR{$VRorig}$N\n" if ($read_name eq "\@6acb2984-5480-415b-96db-da1208fa1a0d");
		while ($init == 0 or $minlen >= $minlen_threshold) {
			next if defined $dontuse->{$VRname};
			$init = 1;
			my $i = -1;
			while (1) {
				my ($VR, $VRrev);
				my $substr = $minlen > length($VRorig) ? length($VRorig) : $minlen;
				if ($i == -1) {
					($VR, $VRrev) = ($VRorig, $VRrevorig);
				}
				else {
					$VR = substr($VRorig, $i, $substr);
					$VRrev = revcomp($VR);
				}
				if ($read_seq =~ /($VR|$VRrev)/) {
					my $VRhash;
					if ($read_seq =~ /($VR)/) {
						$VRhash->{$VRorig} = $VR{$VRorig};
						$VRstrand = "+";
					}
					elsif ($read_seq =~ /($VRrev)/) {
						$VRhash->{$VRrevorig} = $VRname;
						$VRstrand = "-";
					}
					
				   my ($VRfinal0, $VRseqfinal0, $VRstrand20, $VRseqbeg0, $VRseqend0, $beg0, $end0, $totalbc0, $bcarr0, $bcseqarr0, $score0, $bestlength0) = allmatch_check2($read_name, $read_seq, $read_qual, $VRhash);
					$dontuse->{$VRname} = 1;
					if ($score0 > $score) {
					   ($VRfinal, $VRseqfinal, $VRstrand2, $VRseqbeg, $VRseqend, $beg, $end, $totalbc, $bcarr, $bcseqarr, $score, $bestlength) = ($VRfinal0, $VRseqfinal0, $VRstrand20, $VRseqbeg0, $VRseqend0, $beg0, $end0, $totalbc0, $bcarr0, $bcseqarr0, $score0, $bestlength0) = allmatch_check2($read_name, $read_seq, $read_qual, $VRhash) if $score < $score0;
						@bc = @{$bcarr};
						@bcseq = @{$bcseqarr};
						$found = 1 if $score >= $perc_match_threshold;
						$found = 2 if $score >= $perc_match_threshold and $score > 95;
						#if ($read_name eq "\@7954960d-2c5b-4e22-95e5-0ea78e8c452e" and $VRname =~ /VR_VR\-2\-70_Variable_region_2_GCskew/ and $found eq 1) {
						if ($found >= 1) {
							LOG($outLog, "$YW$read_name$N $LCY$VRname$N ${LGN}GOOD$N $LGN$score$N bestlen=$LGN$bestlength$N i=$YW$i to $LGN" . ($i+$substr) . "$N minlen=$LGN$minlen$N\n");
						}
						else {
							$dontuse->{$VRname} = 1;
							LOG($outLog, "$YW$read_name$N $LCY$VRname$N ${LRD}SCORE$N $LGN$score$N bestlen=$LGN$bestlength$N i=$YW$i to $LGN" . ($i+$substr) . "$N minlen=$LGN$minlen$N\n");
						}
					}
						#print "FOUND\n";
					#}
					last if defined $dontuse->{$VRname};
					last if $found eq 1 and not defined $opt_a;
					last if $found eq 2;
				}
				$i += $step_size;
				last if defined $dontuse->{$VRname};
				last if length($VRorig) < $i + $minlen;
				last if $found eq 1 and not defined $opt_a;
				last if $found eq 2;
			}
			$minlen -= $window_subtract_size;
			if (length($VRorig) < $i + $minlen) {
				my $lengthVRorig = length($VRorig);
				#print "last because length VRorig ($lengthVRorig) < pos+minlen=$i+$minlen\n";
			}
			last if length($VRorig) < $i + $minlen;
			last if defined $dontuse->{$VRname};
			last if $found eq 1 and not defined $opt_a;
			last if $found eq 2;
		}
		#print "  while done\n" if ($read_name eq "\@6acb2984-5480-415b-96db-da1208fa1a0d");
		last if $found eq 1 and not defined $opt_a;
		last if $found eq 2;
	}
	#LOG($outLog, "$YW$read_name$N " . scalar(keys %{$dontuse}) . "\n\n");
	if ($VRfinal eq "UNKNOWN") {
		$data{sample}{none} ++;
		$VRseqfinal = join("", ("-") x length($read_seq));
		print $out0 "$read_name\t$beg\t$end\t0\t$VRfinal\t0\t$VRstrand\t$VRseqfinal\t$read_seq\t$VRseqfinal\t$read_seq\n";
		print $out0b "$read_name\t$score\t$bestlength\tUNKNOWN\n";
		LOG($outLog, "$YW$read_name$N $LCY$VRfinal$N ${LRD}UNK$N\n");
		#print "$YW$read_name$N: $LCY$VRfinal$N\n$LGN$read_seq$N\n";
	}
	else {
		$VRfinal = join(";", @bc);
		$VRseqfinal = join(";", @bcseq);
		$data{sample}{$VRfinal} ++;
		$data{samplebc}{$VRseqfinal} ++;
		print $out0b "$read_name\t$score\t$bestlength\t" . join(";", @bc) . "\n";
	}
	#die "good\n";
}
LOG($outLog, "Finished matching on $LGN$linecount$N Reads\n",$opt_q);
close $out0;
close $inFastQ;

$data{total} = 0 if not defined $data{total};

$data{multiplebc} = 0 if not defined $data{multiplebc};
open (my $out1, ">", $VRstatsFile) or die "Cannot write to $VRstatsFile: $!\n";
print $out1 "#total read with more than 1 match: $data{multiplebc}\n";
print $out1 "#match_name\tmatch_perc\tmatch_count\ttotal_read\n";
foreach my $type (sort keys %data) {
	next if $type eq "total";
	next if $type !~ /^(bc1|sample|bc2)$/;
	foreach my $val (sort keys %{$data{$type}}) {
		my $count = $data{$type}{$val};
		my $perc = int($count/$data{total}*1000+0.5)/10;
		print $out1 "$val\t$perc\t$count\t$data{total}\n";
	}
}
close $out1;

LOG($outLog,  "
Output:
$LCY$outFile$N # all reads with VR match
$LCY$outFileshort$N # shorter version of the above
$LCY$logFile$N # log file
$LCY$VRstatsFile$N # stats summary

From $LCY$VRstatsFile$N: 
",$opt_q);

system("cat $VRstatsFile") if not defined $opt_q;

print "\n\n" if not defined $opt_q;
### SUB ROUTINES ###

sub LOG {
   my ($outLog, $text, $STEP, $STEPCOUNT) = @_;
	if (not defined $STEP) {
      print $text;
   }
   if (defined $outLog) {
      print $outLog $text;
   }
}

sub getFilename {
   my ($fh, $type) = @_;
   my $fh0 = $fh;
   $fh =~ s/\/+/\//g;
   my ($first, $last) = (0,0);
   if ($fh =~ /^\//) {
      $first = 1;
      $fh =~ s/^\///;
   }
   if ($fh =~ /\/$/) {
      $last = 1;
      $fh =~ s/\/$//;
   }
   # Split folder and fullname
   my (@splitname) = split("\/+", $fh);
   my $fullname = pop(@splitname);
   my @tempfolder = @splitname;
   my $folder;
   if (@tempfolder == 0) {
      $folder = "./";
   }
   else {
      $folder = join("\/", @tempfolder) . "/";
   }
   #$folder = "./" if $folder eq "/";
#  print "\nFolder = $folder\nFile = $fh, fullname=$fullname\n\n";
   # Split fullname and shortname (dot separated)
        @splitname = split(/\./, $fullname);
        my $shortname = $splitname[0];
   if ($first == 1) {$folder = "/$folder";}
   if ($last == 1) {$fullname = "$fullname/";}
   #print "\nFh=$fh0\nFolder=$folder\nShortname=$shortname\nFullname=$fullname\n\n";
#     print "Retunring shortname\n" if not defined $type;

        return($shortname)          if not defined($type);
        return($folder, $fullname)     if defined($type) and $type =~ /folderfull/;
        return($folder, $shortname)    if defined($type) and $type =~ /folder/;
        return($fullname)        if defined($type) and $type =~ /full/;
        return($folder, $fullname, $shortname)  if defined($type) and $type =~ /all/;
}

sub revcomp {
   my ($sequence) = @_;
   $sequence = uc($sequence);
   $sequence =~ tr/ATGC/TACG/;
   $sequence = reverse($sequence);
   return ($sequence);
}




sub allmatch_check2 {
	my ($read_name, $read_seq, $read_qual, $VRhash) = @_;
	my $VRfinal = "UNKNOWN";
	my $VRseqfinal = "";
	my $VRstrand = ".";
	my ($VRseqbeg, $VRseqend);
	my ($beg, $end) = (-1,-1);
	my $totalbc = 0;
	my @bc;
	my @bcseq;
	my %best;
	$best{score} = -1;
	@{$best{print}} = ();
	my $found = 0;
	#print scalar(keys %{$VRhash}) . "\n";
	foreach my $VRorig (sort keys %{$VRhash}) {
		my $VRname = $VRhash->{$VRorig};
		my $minlen = length($VRorig);
		my $init = 0;
		my $VRrevorig = revcomp($VRorig);
		my $minlen_threshold = length($VRorig); #$minlen_threshold_orig eq 0 ? length($VRorig) : $minlen_threshold_orig;
		if (defined $opt_L) {
			$minlen_threshold = length($VRorig); #int($minlen_threshold_perc / 100 * length($VRorig)+0.5);
		}
		my ($VR, $VRrev) = ($VRorig, $VRrevorig);
		if ($read_seq =~ /($VR|$VRrev)/) {
			if ($found == 0) {
				@{$best{print}} = ();
				delete($best{print});
			}
			my ($VR_seq) = $VR;
			my $seqR = $read_seq;
			if ($read_seq =~ /$VR/) {
				$VR_seq = $VR;
				$seqR = $read_seq;
				($VRseqbeg, $VRseqend) = $read_seq =~ /^(.*)$VR(.*)$/;
				$VRfinal = $VRhash->{$VRorig};
				$VRseqfinal = $VR;
				$VRstrand = "+";
			}
			elsif ($read_seq =~ /$VRrev/) {
				$VR_seq = $VRrev;
				$seqR = $read_seq;
				($VRseqbeg, $VRseqend) = $read_seq =~ /^(.*)$VRrev(.*)$/;
				$VRfinal = $VRhash->{$VRorig};
				$VRseqfinal = $VRrev;
				$VRstrand = "-";
			}
			$beg = length($VRseqbeg);
			$end = length($read_seq) - length($VRseqend);
			my $VRseqbegdot = join("", ("-") x length($VRseqbeg));
			my $VRseqenddot = join("", ("-") x length($VRseqend));
			#if ($totalbc == 1) {
			push(@bc, $VRfinal);
			push(@bcseq, $VRseqfinal);
			$totalbc ++;
			#	$data{sample}{$VRfinal} ++;
			#	$data{samplebc}{$VRseqfinal} ++;
			#}
			if ($totalbc > 1) {
				$data{multiplebc}{$read_name} ++;
			}
			#my $length = $end - $beg;
			my $length = length($VRseqfinal);
			#print "$LCY$VRfinal$N $LPR$VRstrand$N $LGN$best{score} $length/$length$N\n";
			$best{score} = perc($length/length($VRorig));
			$best{len} = $length;
			push(@{$best{print}}, "$read_name\t$beg\t$end\t$best{score}\t$VRfinal\t$totalbc;0;$length\t$VRstrand\t$VRseqbegdot$VRseqfinal$VRseqenddot\t$seqR\t$VR_seq\t$read_seq\t$read_qual");
			#print "$YW$read_name$N: $LCY$VRfinal$N ($beg-$end $VRstrand)\n$LGN$read_seq$N\n$LPR$VRseqbegdot$LGN$VRseqfinal$N$LPR$VRseqenddot$N\n\n";
			#print "   i=$LGN$i$N: LCY$VRfinal$N $LGN$VRhash->{$VRorig}$N\n" if ($read_name eq "\@6acb2984-5480-415b-96db-da1208fa1a0d");
			$found = 1;
		}
		else {
			my @strand = ("+");#, "-");
			my @VRs = ($VR);#, $VRrev);
			for (my $j = 0; $j < @VRs; $j++) {
				my $VR_seq = $VRs[$j];

				my $input = ">$read_name\n$read_seq\n>$VRname\n$VR_seq";
				my $outFile = ".VRcheckertemp/$fqFilename\_$VRFilename\_l$minlen_threshold_orig.temp" if defined $opt_l;
				   $outFile = ".VRcheckertemp/$fqFilename\_$VRFilename\_L$minlen_threshold_orig.temp" if defined $opt_L;
				open (my $out, ">", $outFile) or die;
				print $out "$input\n";
				close $out;
				my $cmd = "cat $outFile | clustalo-1.2.4-Ubuntu-x86_64 -i -";
				#print "$cmd\n";
				my $cmdout = `$cmd`;
				my ($defR, $seqR, $defV, $seqV) = parse_clustalo_res($cmdout);
				my ($matfinal,$mat,$lenfinal,$begfinal,$endfinal,$seqfinal) = count_match($defR, $seqR, $defV, $seqV, $minlen_threshold);
				#print "$LCY$VRname$N $LPR$strand[$j]$N $LGN$matfinal $mat/$lenfinal$N\n";
				if ($best{score} <= $matfinal) {
					if ($best{score} < $matfinal) {
						$best{totalbc} = 0;
						$best{res} = ();
						$best{print} = ();
						$best{bc} = ();
						$best{bcseq} = ();
						$best{len} = 0;
						delete($best{res});
						delete($best{print});
						delete($best{bc});
						delete($best{bcseq});
					}
					$best{score} = $matfinal;
					$VRfinal = $VRhash->{$VRorig};
					$VRseqfinal = $seqfinal;
					$VRstrand = $strand[$j];
					($VRseqbeg) = $seqV =~ /^(.{$begfinal})/; 
					($VRseqend) = $seqV =~ /^.{$endfinal}(.+)$/; 
					#if (not defined $VRseqbeg) {
					#	print "$LGN$read_name $VRname $LCY$seqV$N\nbegfinal=$begfinal\n";
					#}
					#if (not defined $VRseqend) {
					#	print "$LGN$read_name $VRname $LCY$seqV$N\nendfinal=$endfinal\n";
					#}
					$VRseqbeg = "" if not defined $VRseqbeg;
					$VRseqend = "" if not defined $VRseqend;
					my $VRseqbegdot = $VRseqbeg;
					my $VRseqenddot = $VRseqend;
					$beg = $begfinal;
					$end = $endfinal;
					my $length = $lenfinal;
					$best{len} = $length;
					$best{totalbc} ++;
					push(@{$best{bc}}, $VRfinal);
					push(@{$best{bcseq}}, $VRseqfinal);
					push(@{$best{res}}, "$VRfinal\t$VRseqfinal\t$VRstrand\t$VRseqbeg\t$VRseqend\t$beg\t$end\t$VRseqbegdot");
					push(@{$best{print}}, "$read_name\t$beg\t$end\t$best{score}\t$VRfinal\tTOTALBC;1;$length\t$VRstrand\t$seqR\t$YW$VRseqbegdot$LCY$VRseqfinal$YW$VRseqenddot$N\t$seqR\t$seqV\t$read_seq\t$VR_seq\t$read_qual");
				}
				#print "$LCY$res$N\n";
			#print $out0 "$read_name\t$beg\t$end\t$VRfinal\t$totalbc;$i;$length\t$VRstrand\t$VRseqbegdot$VRseqfinal$VRseqenddot\t$seqR\t$VR_seq\t$read_seq\t$read_qual\n";
			}
		}
	}
	if ($found ne 1) {
		my $totalbc = $best{totalbc};
		@bc = @{$best{bc}};
		@bcseq = @{$best{bcseq}};
		($VRfinal, $VRseqfinal, $VRstrand, $VRseqbeg, $VRseqend, $beg, $end) = split("\t", $best{res}[0]);
		if ($totalbc > 1) {
			$data{multiplebc}{$read_name} ++;
		}
		for (my $i = 0; $i < @{$best{print}}; $i++) {
			$best{print}[$i] =~ s/\tTOTALBC;/\t$totalbc;/;
			print $out0 "$best{print}[$i]\n";
		}
	}
	else {
		for (my $i = 0; $i < @{$best{print}}; $i++) {
			print $out0 "$best{print}[$i]\n";
		}
	}
	return($VRfinal, $VRseqfinal, $VRstrand, $VRseqbeg, $VRseqend, $beg, $end, $totalbc, \@bc, \@bcseq, $best{score},$best{len});
}
sub parse_clustalo_res {
	my ($resline) = @_;
	my $res;
	my @line = split("\n", $resline);
	my ($def, $seq);
	my $seqcount = 0;
	my @res;
	for (my $i = 0; $i < @line; $i++) {
		chomp($line[$i]);
		if ($line[$i] =~ /^>/) {
			if (defined $seq) {
				#print "$YW$def\n$seq$N\n";
				@res = (@res, $def, $seq);
			}
			$seq = "";
			$def = $line[$i];
		}
		else {
			$seq .= $line[$i];
		}
	}
	if (defined $seq) {
		#print "$YW$def\n$seq$N\n";
		@res = (@res, $def, $seq);
	}
	return(@res);
}

sub count_match {
	my ($defR, $seqsR, $defV, $seqsV, $minlen_threshold) = @_;
	#print "$LCY$defR$N\n$seqsR\n$LCY$defV$N\n$seqsV\n";
	my @seqsR = split("", $seqsR);
	my @seqsV = split("", $seqsV);
	my ($edgebeg) = $seqsV =~ /^([\-]*)[A-Z]/;
	$edgebeg = defined $edgebeg ? length($edgebeg) : 0;
	my ($edgeend) = $seqsV =~ /^(.+[A-Z])[\-]*$/;
	$edgeend = defined $edgeend ? length($edgeend) : 0;
	my @res = ("U") x @seqsV;
	for (my $i = 0; $i < $edgebeg; $i++) {
		$res[$i] = " ";
	}
	for (my $i = $edgeend; $i < @seqsV; $i++) {
		$res[$i] = " ";
	}
	#print "edge=$LGN$edgebeg-$edgeend$N\n";
	my @mat;
	my $mat = 0;
	my $mis = 0;
	my $total = 0;
	my $last = 0;
	my @match;
	my $seqsVthres = $seqsV;
	$seqsVthres =~ s/\-//g;
	$seqsVthres = length($seqsVthres);
	#my $minlen_threshold = $minlen_threshold_orig eq 0 ? ($seqsVthres) : $minlen_threshold_orig;
	#if (defined $opt_L) {
	#	$minlen_threshold = int($minlen_threshold_perc / 100 * ($seqsVthres)+0.5);
	#}
	for (my $i = 0; $i < @seqsV; $i++) {
		my $nucR = $seqsR[$i]; # read
		my $nucV = $seqsV[$i]; # VR
		# AAA---TTT---AAA
		# --ATGCTTTGGG-A-
		# EEMIIIMMMIIIDME
		# E = 
		# D = read N VR -
		# I = read - VR N
		# 
		
		if ($res[$i] eq " ") {
			$mat[$i] = 0;
		}
		else {
			$total ++;
			if ($nucR eq "-" and $nucV eq "-") {
				print "$LCY$defR$N\n$seqsR\n$LCY$defV$N\n$seqsV\n";
				die "Fatal error: both nucRead and nucVR are '-'\n";
				
			}
			elsif ($nucR eq "-") {
				$res[$i] = "i";
				$mat[$i] = 0;
			}
			elsif ($nucV eq "-") {
				$res[$i] = ".";
				$mat[$i] = 0;
			}
			elsif ($nucR eq $nucV) {
				$res[$i] = "M";
				$mat[$i] = 1;
				$mat ++;
			}
			elsif ($nucR ne $nucV) {
				$res[$i] = "m";
				$mat[$i] = 0;
			}
			else {
				die "Fatal error: something unaccounted (nucR=$LCY$nucR$N, nucV=$LGN$nucV$N)\n";
			}
		}
	}
	my ($matfinal, $begfinal, $endfinal, $lenfinal, $seqfinal)=(0,0,0);
	my ($matfinal2, $begfinal2, $endfinal2, $lenfinal2, $seqfinal2)=(0,0,0);

	$matfinal = $total == 0 ? 0 : perc($mat/$total);
	$begfinal = $edgebeg;
	$endfinal = $edgeend;
	$lenfinal = ($endfinal-$begfinal);
	($seqfinal) = $seqsV =~ /^.{$begfinal}(.{$lenfinal}).*$/;
	my $mat2 = 0;
	for (my $i = 0; $i < @mat-$minlen_threshold; $i++) {
		my $currbeg = $i;
		my $currend = $i + $minlen_threshold - 1;
		#print "$currbeg-$currend\n" if $i == 0;
		my @currmat = @mat[$currbeg..$currend];
		my $currmat = perc(sum(\@currmat) / $minlen_threshold);
		if ($currmat > $matfinal2) {
			$mat2 = sum(\@currmat);
			$matfinal2 = $currmat;
			$begfinal2 = $i;
			$endfinal2 = $i+$minlen_threshold;
			$lenfinal2 = $endfinal2-$begfinal2;
			($seqfinal2) = $seqsV =~ /^.{$begfinal2}(.{$lenfinal2}).*$/;
		}
	}
	my $matfinal1 = $matfinal;
	my $begfinal1 = $begfinal;
	my $endfinal1 = $endfinal;
	my $lenfinal1 = $lenfinal;
	my $seqfinal1 = $seqfinal;
	my $mat1 = $mat;
	if ($matfinal2 > $matfinal) {
		$matfinal = $matfinal2;
		$begfinal = $begfinal2;
		$endfinal = $endfinal2;
		$lenfinal = $lenfinal2;
		$seqfinal = $seqfinal2;
		$mat = $mat2;
	}
	#print "$seqsR\n$seqsV\n";
	#print join("", @res) . "\n";
	#print "mat=$mat, total=$total, minlen_threshold=$minlen_threshold, seqsvthres=$seqsVthres\n";
	#print "match final1=$LGN$matfinal1\%$N ($LGN$mat1$N/$LGN$lenfinal1$N) ($LCY$begfinal1-$endfinal1$N)\n";
	#print "match final2=$LGN$matfinal2\%$N ($LGN$mat2$N/$LGN$lenfinal2$N) ($LCY$begfinal2-$endfinal2$N) minlen_threshold=$LGN$minlen_threshold$N\n";
	#print "match final=$LGN$matfinal\%$N ($LGN$mat$N/$LGN$lenfinal$N) len=$LGN$lenfinal$N ($LCY$begfinal-$endfinal$N)\n";
	return($matfinal,$mat,$lenfinal,$begfinal,$endfinal,$seqfinal);
}
sub perc {
	my ($num) = @_;
	return(int(1000*$num+0.5)/10);
}
sub mean {
	my ($array) = @_;
	my $mean = 0;
	return 0 if (@{$array} == 0);
	return 0 if $array !~ /ARRAY/;
	for (my $i = 0; $i < @{$array}; $i++) {
		$mean += $array->[$i] / @{$array};
	}
	return($mean);
}

sub sum {
	my ($array) = @_;
	my $sum = 0;
	return 0 if (@{$array} == 0);
	return 0 if $array !~ /ARRAY/;
	for (my $i = 0; $i < @{$array}; $i++) {
		$sum += $array->[$i];
	}
	return($sum);
}

