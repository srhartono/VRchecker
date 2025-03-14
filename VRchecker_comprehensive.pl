#!/usr/bin/perl

use strict; use warnings; use Getopt::Std; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_i $opt_b $opt_s $opt_w $opt_l $opt_o $opt_H $opt_q $opt_L $opt_t $opt_a $opt_c $opt_0 $opt_A $opt_T);
getopts("vi:b:s:w:l:o:HqL:t:aAc:0T:");

my $libPath;
BEGIN {
	$libPath = dirname(dirname abs_path $0) . '/VRchecker/';
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

my $default_step_size             = 20;
my $default_window_subtract_size  = 20;
my $default_minlen_threshold_orig = 20;
my $default_minlen_threshold_perc = 10;
my $default_perc_match_threshold  = 75;
my $default_outDir                = "./";
my $default_clustaloFile          = "$libPath/clustalo-1.2.4-Ubuntu-x86_64";

my ($fqFile, $VRFile) = ($opt_i, $opt_b);
my $usage = "

Usage: $YW$0$N -i $LCY<fqFile.fastq.gz>$N -b $LPR<VRFile.fa>$N

-T: length for checking match threshold [VRlength]
-t: percent match threshold [75]
-o: output dir [currrent directory]
-q: quiet, don't print anything

-c: clustalo aligner
    [default: $LGN$libPath$N/${YW}clustalo-1.2.4-Ubuntu-x86_64$N]

Seeding match parameters:
-s: step size [20]
-w: window subtract size [20]
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
-w: window subtract size [20]
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

die $usage_long if defined $opt_H;
die $usage unless defined $opt_i and defined $opt_b;
die "fastq file -i $LCY$fqFile$N doesn't exist!\n" if not -e $fqFile;
die "query file -b $LCY$VRFile$N doesn't exist!\n" if not -e $VRFile;

my $step_size             = defined $opt_s ? $opt_s : $default_step_size;
my $window_subtract_size  = defined $opt_w ? $opt_w : $default_window_subtract_size;
my $minlen_threshold_orig = defined $opt_l ? $opt_l : $default_minlen_threshold_orig;
my $minlen_threshold_perc = defined $opt_L ? $opt_L : $default_minlen_threshold_perc;
my $PERC_MATCH_THRESHOLD  = defined $opt_t ? $opt_t : $default_perc_match_threshold;
my $outDir                = defined $opt_o ? $opt_o : $default_outDir;
my $clustaloFile          = defined $opt_c ? $opt_c : $default_clustaloFile;
my $quiet = $opt_q;

die "\n\n${LRD}ERROR$N: Cannot find clustalo file $YW$clustaloFile$N anywhere!\nTry putting $libPath/clustalo-1.2.4-Ubuntu-x86_64 to \$PATH\n\n" if not -e $clustaloFile;

my %VR;
system("mkdir -p $outDir") if not -d $outDir;
$outDir =~ s/\/+$//;
my ($folder1, $fqFilename) = getFilename($fqFile, "folderfull");
my ($folder2, $VRFilename) = getFilename($VRFile, "folderfull");
my $logFile = "$outDir/$fqFilename\_$VRFilename\_l$minlen_threshold_orig\_w$window_subtract_size\_s$step_size.comprehensive.log";
my $outFile = "$outDir/$fqFilename\_$VRFilename\_l$minlen_threshold_orig\_w$window_subtract_size\_s$step_size.comprehensive.tsv";
my $outFileshort = "$outDir/$fqFilename\_$VRFilename\_l$minlen_threshold_orig\_w$window_subtract_size\_s$step_size.comprehensive.short.tsv";
my $VRstatsFile = "$outDir/$fqFilename\_$VRFilename\_l$minlen_threshold_orig\_w$window_subtract_size\_s$step_size.comprehensive.stats";
system("mkdir -p $outDir/.VRchecker_comprehensive/") if not -d "$outDir/.VRchecker_comprehensive/";
my $clustaloutFile = "$outDir/.VRchecker_comprehensive/$fqFilename\_$VRFilename\_t$PERC_MATCH_THRESHOLD\_l$minlen_threshold_orig\_L$minlen_threshold_perc\_w$window_subtract_size\_s$step_size.comprehensive.clustalo.temp";
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
-t perc match threshold  : $LGN$PERC_MATCH_THRESHOLD$N

Added VRchecker directory $LCY$libPath$N to perl ${LCY}\@INC$N
Clustal output temp file : $LCY$clustaloutFile$N

",$opt_q);
LOG($outLog, "\n#VR sequences:\n",$opt_q);

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
if (defined $opt_0) { #debug
	LOG($outLog, "\n" . date() . "${LGN}Debug Successful!!$N\n\n");
	exit 0;
}

LOG($outLog, "\n# matching query sequence to fastq sequence\n",$opt_q);
my $linecount = 0;
my %data = ("total" => 0, "multiplebc" => 0);
#$data{multiplebc} = 0 if not defined $data{multiplebc};
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
print $out0 "#read_name\tbeg\tend\tbestscore\tbestVRstrand\tVRname\tTOTALBC;FOUND;bestlength\tVRseq\tseqR\tVRseqbegbestVRseqVRseqend\tseqR\tseqV\tread_seq\tread_qual\n";
#read_name\tmatch_start\tmatch_end\tmatch_name\tmatch_number\tmatch_strand\tmatch_seq\tread_seq\tread_qual\n";
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
	seed_checker($read_name, $read_seq, $read_qual);
}
LOG($outLog, "Finished matching on $LGN$linecount$N Reads\n",$opt_q);
close $out0;
close $inFastQ;

my $VRstats;
$VRstats .= "##total read with more than 1 match: $data{multiplebc}\n";
$VRstats .= "#match_name\tmatch_perc\tmatch_count\ttotal_read\n";
foreach my $type (sort keys %data) {
	next if $type eq "total";
	next if $type !~ /^(bc1|sample|bc2)$/;
	foreach my $val (sort keys %{$data{$type}}) {
		my $count = $data{$type}{$val};
		my $perc = int($count/$data{total}*1000+0.5)/10;
		$VRstats .= "$val\t$perc\t$count\t$data{total}\n";
	}
}

open (my $out1, ">", $VRstatsFile) or die "Cannot write to $VRstatsFile: $!\n";
print $out1 $VRstats;
close $out1;
my $md5sum = `md5sum $VRstatsFile`;
LOG($outLog,  "
Output:
$LCY$outFile$N # all reads with VR match
$LCY$outFileshort$N # shorter version of the above
$LCY$logFile$N # log file
$LCY$VRstatsFile$N # stats summary

From $LCY$VRstatsFile$N: 
$VRstats
\n\n",$opt_q);


LOG($outLog, "$LGN$md5sum$N\n\n");
### SUB ROUTINES ###

sub init_temphash {
	my $temp;
	$temp->{read_name} = "";
	$temp->{read_seq} = "";
	$temp->{read_qual} = "";
	$temp->{VR} = "UNKNOWN";
	$temp->{VRseq} = "UNKNOWN";
	$temp->{bestVRname} = "UNKNOWN";
	$temp->{bestVRseq} = "UNKNOWN";
	$temp->{bestVRstrand} = "UNKNOWN";
	$temp->{VRlength} = "";
	$temp->{VRstrand} = ".";
	$temp->{VRseqbeg} = "";
	$temp->{VRseqend} = "";
	$temp->{bestbeg} = -1;
	$temp->{bestend} = -1;
	$temp->{beg} = -1;
	$temp->{end} = -1;
	$temp->{total} ++;
	@{$temp->{bcname}} = ();
	@{$temp->{bcseq}} = ();
	@{$temp->{bcstrand}} = ();
	$temp->{bestscore} = 0;
	$temp->{bestlength} = 0;
	%{$temp->{VRchecked}} = ();
	return($temp);
}
sub arrcopy {
	my ($arr1) = @_;
	my $arr2;
	for (my $i = 0; $i < @{$arr1}; $i++) {
		if ($arr1->[$i] =~ /ARRAY/) {
			$arr2->[$i] = arrcopy($arr1->[$i])
		}
		elsif ($arr1->[$i] =~ /HASH/) {
			$arr2->[$i] = hashcopy($arr1->[$i])
		}
		else {
			$arr2->[$i] = $arr1->[$i];
		}
	}
	return($arr2);
}
sub hashcopy {
	my ($hash1) = @_;
	my $hash2;
	foreach my $key (keys %{$hash1}) {
		if ($hash1->{$key} =~ /ARRAY/) {
			$hash2->{$key} = arrcopy($hash1->{$key})
		}
		elsif ($hash1->{$key} =~ /HASH/) {
			$hash2->{$key} = hashcopy($hash1->{$key})
		}
		else {
			$hash2->{$key} = $hash1->{$key};
		}
	}
	return($hash2);
}
sub seed_checker {
	my ($read_name, $read_seq, $read_qual) = @_;
	my $temp = init_temphash();
	$temp->{read_name} = $read_name;
	$temp->{read_seq} = $read_seq;
	$temp->{read_qual} = $read_qual;
	my $last = 0;
	$data{total} ++;
	#print "\n$YW$read_name$N\n";
	foreach my $VRseq (sort keys %VR) {		
		last if $last eq 2;
		$temp->{VRseq} = $VRseq;
		#LOG($outLog, " $LCY$VR{$VRseq}$N\n");
		$temp->{VRname} = $VR{$VRseq};
		$temp->{VRlength} = length($VRseq);
		die if defined $temp->{VRchecked}{$temp->{VRname}};
		$temp->{found} = 0;
		$temp->{minlen_threshold} = defined $opt_L ? int($temp->{minlen_threshold_perc} / 100 * length($VRseq)+0.5) : $minlen_threshold_orig;
	 	$temp->{minlen_threshold} = min(length($VRseq),max(0,$temp->{minlen_threshold}));
		for ($temp->{minlen} = $temp->{VRlength}; $temp->{minlen} >= $temp->{minlen_threshold}; $temp->{minlen} -= $window_subtract_size) {
			last if $last eq 1;
			#LOG($outLog, "   $temp->{minlen}\n");
		#while ($last eq 0 and ($init == 0 or $minlen >= $temp->{minlen_threshold})) {
			#$init = 1;
			for (my $posbeg = 0; $posbeg <= $temp->{VRlength} - $temp->{minlen}; $posbeg += $step_size) {
				last if $last eq 1;
				my $posend = $posbeg + $temp->{minlen};
				#LOG($outLog, "     $posbeg-$posend\n");
				$temp->{VRchunk} = substr($temp->{VRseq}, $posbeg, $temp->{minlen});
				#print "substr " . length($temp->{VRseq}) . ", $posbeg, $temp->{minlen}\n";
				$temp->{VRstrand} = "+";
				if (revcomp($read_seq) !~ /$temp->{VRchunk}/i and $read_seq !~ /$temp->{VRchunk}/i) {
					next;
				}
				#elsif (revcomp($read_seq) =~ /$temp->{VRchunk}/i) {
				#	$temp->{VRstrand} = "-";
				#	$temp->{VRseq} = revcomp($temp->{VRseq});
				#}

				for (my $k = 0; $k < 2; $k++) {
					$temp->{VRseq} = revcomp($temp->{VRseq}) if $k == 1;
					$temp->{VRstrand} = "-" if $k == 1;
					my $prev_bestscore = $temp->{bestscore};
					$temp = allmatch_check2($temp);
					#LOG($outLog, "best=$temp->{bestscore} $temp->{bestVRname} beststrand=$temp->{bestVRstrand}\n");
					if ($prev_bestscore < $temp->{bestscore}) {
						$temp->{found} = 1 if $temp->{bestscore} >= $PERC_MATCH_THRESHOLD;
						$temp->{found} = 2 if $temp->{bestscore} >= $PERC_MATCH_THRESHOLD	 and $temp->{bestscore} > 95;
						if ($temp->{found} >= 1) {
							#LOG($outLog, "- $YW$read_name$N $LCY$temp->{VRname}$N ${LGN}GOOD$N $LGN$temp->{bestscore}$N $LGN$temp->{bestVRstrand}$N bestlen=$LGN$temp->{bestlength}$N i=$LGN$posbeg-$posend$N ($LGN$temp->{minlen}$N bp) minlen=$LGN$temp->{minlen}$N\n",$opt_q);
						}
						else {
							#LOG($outLog, "- $YW$read_name$N $LCY$temp->{VRname}$N ${LRD}LOWSCORE$N $LGN$temp->{bestscore}$N $LRD$temp->{bestVRstrand}$N bestlen=$LGN$temp->{bestlength}$N i=$LGN$posbeg-$posend$N ($LGN$temp->{minlen}$N bp) minlen=$LGN$temp->{minlen}$N\n",$opt_q);
						}
						$temp->{VRchecked}{$temp->{VRname}} = 1;
					}
					else {
						#LOG($outLog, "- $YW$read_name$N $LCY$temp->{VRname}$N ${LGN}KEPT$N $LGN$temp->{bestscore}$N $LGN$temp->{bestVRstrand}$N bestlen=$LGN$temp->{bestlength}$N i=$LGN$posbeg-$posend$N ($LGN$temp->{minlen}$N bp) minlen=$LGN$temp->{minlen}$N\n",$opt_q);
					}
					$last = 2 if $temp->{found} eq 1 and not defined $opt_a and not defined $opt_A;
					$last = 2 if $temp->{found} eq 2 and not defined $opt_A;
				}
				$temp->{VRchecked}{$temp->{VRname}} = 1;
				$last = 1 if defined $temp->{VRchecked}{$temp->{VRname}};
				#$last = 1 if length($VRseq) < $posbeg + $temp->{minlen};
			}
			#$temp->{minlen} -= $window_subtract_size;
			#$last = 1 if length($VRseq) < $posbeg + $temp->{minlen};
			$last = 2 if $temp->{found} eq 1 and not defined $opt_a and not defined $opt_A;
			$last = 2 if $temp->{found} eq 2 and not defined $opt_A;
			$last = 1 if defined $temp->{VRchecked}{$temp->{VRname}};
		}
		$last = 2 if $temp->{found} eq 1 and not defined $opt_a and not defined $opt_A;
		$last = 2 if $temp->{found} eq 2 and not defined $opt_A;
	}
	if ($temp->{found} == 0) {#bestVRname} eq "UNKNOWN") {
		$data{sample}{none} ++;
		if ($temp->{bestVRname} ne "UNKNOWN") {
			my $totalbc = @{$temp->{bcname}};
			for (my $i = 0; $i < @{$temp->{bestprint}}; $i++) {
				$temp->{bestprint}[$i] =~ s/\tTOTALBC;/\t$totalbc;/;
				$temp->{bestprint}[$i] =~ s/;FOUND;/;$temp->{found};/;
				print $out0 "$temp->{bestprint}[$i]\n";
			}
		}
		else {
			print $out0 "$temp->{read_name}\t$temp->{beg}\t$temp->{end}\t$temp->{bestscore}\t$temp->{bestVRstrand}\t$temp->{bestVRname}\t0;0;$temp->{bestlength}\t$temp->{bestVRseq}\t$temp->{read_seq}\t$temp->{VRseqbeg}$LCY$temp->{bestVRseq}$temp->{VRseqend}$N\t$temp->{read_seq}\t$temp->{bestVRseq}\t$temp->{read_seq}\t$temp->{read_qual}\n";
		}
		print $out0b "$read_name\t$temp->{bestscore}\t$temp->{bestlength}\t$temp->{bestVRname}\t$temp->{bestVRstrand}\n";
		LOG($outLog, "$YW$read_name$N $LGN$temp->{bestscore}$N $LGN$temp->{bestbeg}-$temp->{bestend} ($temp->{bestVRstrand})$N $LGN$temp->{bestlength}$N $LCY$temp->{bestVRname}$N\n",$opt_q);
	}
	else {
		$temp->{bestVRname} = join(";", @{$temp->{bcname}});
		$temp->{bestVRseq} = join(";", @{$temp->{bcseq}});
		$temp->{bestVRstrand} = join(";", @{$temp->{bcstrand}});
		$temp->{bestVRmat} = join(";", @{$temp->{bcmat}});
		$data{sample}{$temp->{bestVRname}} ++;
		$data{samplebc}{$temp->{bestVRseq}} ++;
		$data{samplestrand}{$temp->{bestVRstrand}} ++;
		$data{samplemat}{$temp->{bestVRmat}} ++;
		my $totalbc = @{$temp->{bcname}};

		for (my $i = 0; $i < @{$temp->{bestprint}}; $i++) {
			$temp->{bestprint}[$i] =~ s/\tTOTALBC;/\t$totalbc;/;
			$temp->{bestprint}[$i] =~ s/;FOUND;/;$temp->{found};/;
			print $out0 "$temp->{bestprint}[$i]\n";
		}
		print $out0b "$read_name\t$temp->{bestscore}\t$temp->{bestlength}\t$temp->{bestVRname}\t$temp->{bestVRstrand}\n";
		LOG($outLog, "$YW$read_name$N $LGN$temp->{bestscore}$N $LGN$temp->{bestbeg}-$temp->{bestend} ($temp->{bestVRstrand})$N $LGN$temp->{bestlength}$N $LCY$temp->{bestVRname}$N\n",$opt_q);
		if ($totalbc > 1) {
			$data{multiplebc}{$temp->{read_name}} ++;
		}
	}
	#foreach my $type (sort keys %data) [
	#	foreach my $VRnamefinal (sort keys %{$data{sample}}) {
	#		print "$VRnamefinal: $data{sample}{$VRnamefinal}\n";
	#	}
	#}
}

sub allmatch_check2 {
	my ($temp) = @_;
	my ($defR, $seqR, $defV, $seqV) = run_clustalo($temp);
	run_clustalo($temp,1) if not defined $seqR;
	#print "$seqR\n$seqV\n";
	die "Allmatch check2 run_clostalo seqR undef\n" if not defined $seqR;
	my ($matfinal,$mat,$lenfinal,$begfinal,$endfinal,$seqfinal) = count_match($defR, $seqR, $defV, $seqV, $temp->{VRlength}, $temp->{minlen_threshold});
	#LOG($outLog, "$temp->{VRname} $temp->{VRstrand}: $matfinal\n");
	if ($temp->{bestscore} <= $matfinal) {
		return($temp) if $temp->{bestscore} == $matfinal and $temp->{VRname} eq $temp->{bestVRname};
		if ($temp->{bestscore} < $matfinal) {
			$temp->{bestscore} = $matfinal;
			@{$temp->{bestprint}} = ();
			@{$temp->{bcname}} = ();
			@{$temp->{bcseq}} = ();
			@{$temp->{bcstrand}} = ();
			@{$temp->{bcmat}} = ();
			@{$temp->{bclength}} = ();
		}
		$temp->{bestVRname} = $temp->{VRname};
		$temp->{bestVRseq} = $seqfinal;
		$temp->{bestVRstrand} = $temp->{VRstrand};
		$temp->{bestVRmat} = $mat;
		($temp->{beg}, $temp->{end}) = ($begfinal, $endfinal);
		($temp->{bestbeg}, $temp->{bestend}) = ($begfinal, $endfinal);
		($temp->{VRseqbeg}) = $begfinal > 0 ? $seqV =~ /^(.{$begfinal})/ : "";
		($temp->{VRseqend}) = $endfinal < length($seqV) ? $seqV =~ /^.{$endfinal}(.*)$/ : "";
		$temp->{bestlength} = $lenfinal;
		push(@{$temp->{bcname}}, $temp->{bestVRname});
		push(@{$temp->{bcseq}}, $temp->{bestVRseq});
		push(@{$temp->{bcstrand}}, $temp->{bestVRstrand});
		push(@{$temp->{bcmat}}, $temp->{bestVRmat});
		push(@{$temp->{bclength}}, $temp->{bestlength});
		push(@{$temp->{bestprint}}, "$temp->{read_name}\t$temp->{beg}\t$temp->{end}\t$temp->{bestscore}\t$temp->{bestVRstrand}\t$temp->{VRname}\tTOTALBC;FOUND;$temp->{bestlength}\t$temp->{VRseq}\t$seqR\t$LPR$temp->{VRseqbeg}$LCY$temp->{bestVRseq}$LGN$temp->{VRseqend}$N\t$seqR\t$seqV\t$temp->{read_seq}\t$temp->{read_qual}");
	}
	return($temp);#best{score}, $best{length}, $VRfinal, $VRseqfinal, $VRseqbeg, $VRseqend, $beg, $end, $totalbc, \@bc, \@bcseq);
}
sub run_clustalo {
	my ($temp, $debug) = @_;
	my $read_name_edited = $temp->{read_name};
		$read_name_edited =~ s/^\@//;
	my $input = ">$read_name_edited\n$temp->{read_seq}\n>$temp->{VRname}\n$temp->{VRseq}";
	open (my $out, ">", $clustaloutFile) or die "allmatch_check2:: Failed to write to outFile $clustaloutFile: $!\n\n";
	print $out "$input\n";
	close $out;
	print "input=$input\n" if defined $debug;
	return if defined $debug;
	my $cmdout = `cat $clustaloutFile | $clustaloFile -i -`;
	return(parse_clustalo($cmdout));
}
sub parse_clustalo {
	my ($cmdout) = @_;
	my $res;
	my @line = split("\n", $cmdout);
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
	my ($defR, $seqsR, $defV, $seqsV, $VRlength, $minlen_threshold) = @_;
	my $denom_length = defined $opt_T ? $opt_T : $VRlength;
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
	for (my $i = 0; $i < @mat-$denom_length; $i++) {
		my $currbeg = $i;
		my $currend = $i + $denom_length - 1;
		#print "$currbeg-$currend\n" if $i == 0;
		my @currmat = @mat[$currbeg..$currend];
		my $currmat = perc(sum(\@currmat) / $denom_length);
		if ($currmat > $matfinal2) {
			$mat2 = sum(\@currmat);
			$matfinal2 = $currmat;
			$begfinal2 = $i;
			$endfinal2 = $i+$denom_length;
			$lenfinal2 = $endfinal2-$begfinal2;
			($seqfinal2) = $seqsV =~ /^.{$begfinal2}(.{$lenfinal2}).*$/;
			$seqfinal2 =~ s/([A-Z])[\-]*$/$1/;
			$lenfinal2 = length($seqfinal2);
			$endfinal2 = $begfinal2 + $lenfinal2;
			$seqfinal2 =~ s/^[\-]*([A-Z])/$1/;
			$lenfinal2 = length($seqfinal2);
			$begfinal2 = $endfinal2 - $lenfinal2;
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
	#print "mat=$mat, total=$total, $denom_length=$denom_length, seqsvthres=$seqsVthres\n";
	#print "match final1=$LGN$matfinal1\%$N ($LGN$mat1$N/$LGN$lenfinal1$N) ($LCY$begfinal1-$endfinal1$N)\n";
	#print "match final2=$LGN$matfinal2\%$N ($LGN$mat2$N/$LGN$lenfinal2$N) ($LCY$begfinal2-$endfinal2$N) $denom_length=$LGN$denom_length$N\n";
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

sub date {
	my ($add, $color) = @_;
	($color = $add and undef $add) if defined $add and $add =~ /^(color=|color|col=|col|y|c)/i;
	$add = 0 if not defined $add;
	return ("[$YW" . getDate($add) . "$N]: ") if not defined $color;
	return ("[" . getDate($add) . "]: ") if defined $color;
}

sub getDate {
	my ($add) = @_;
	$add = 0 if not defined $add;
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time + $add); $year += 1900;
	my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
	my $date = "$year $mday $months[$mon] $hour:$min:$sec";
	return($date);
}

sub max {
	my (@arr) = @_;
	die "array broken:\n" . join("\n", @arr) . "\nEND\n\n" if not defined $arr[1];
	@arr = sort {$b <=> $a} @arr;
	return($arr[0]);
}
sub min {
	my (@arr) = @_;
	die "array broken:\n" . join("\n", @arr) . "\nEND\n\n" if not defined $arr[1];
	@arr = sort {$a <=> $b} @arr;
	return($arr[0]);
}


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

