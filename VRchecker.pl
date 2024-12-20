#!/usr/bin/perl

use strict; use warnings; use Getopt::Std; use Cwd qw(abs_path); use File::Basename qw(dirname);
use vars qw($opt_v $opt_i $opt_b $opt_s $opt_w $opt_l $opt_o $opt_H $opt_q);
getopts("vi:b:s:w:l:o:Hq");

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
Usage: $YW$0$N -i $LCY<fqFile>$N -b $LPR<VRFile>$N

-o: output dir [currrent directory]
-l: min match length threshold [100]

Examples:

- Exact match only (assuming VR sequence length is 200bp)
$YW$0$N -i ${LCY}001_pFC8-VR-TAC-2x-LAC-Mini-1_reads.fastq.gz$N -b ${LGN}VR.fa$N -l 200

- At least 100bp of VR sequence match:
$YW$0$N -i ${LCY}001_pFC8-VR-TAC-2x-LAC-Mini-1_reads.fastq.gz$N -b ${LGN}VR.fa$N -l 100

More options & explanation do -H

";

my $usage_long = "

Options for quick dirty sliding-window based non-exact matching (blast takes too long)
-s: step size [20]
-w: window subtract size [25]


If the whole VR length doesn't match (200bp) it'll try to reduce the window size by [-w] 
and do sliding window matching with step size of [-s]
It'll stop when finding a match or when window size is less than min length threshold [-l, default 100]

Example: -s step size 20, -w window 25, min match length threshold$LGN 100$N

Initial VR sequence size is 200, try match to the read sequence, if not match then:
window = 200 - 25 = 175 bp -> substring position 0-175 and 20-195 of VR sequence, then try match to read sequence, if not:
window = 175 - 25 = 150 bp -> substring position 0-150, 20-170, and 40-190 of VR sequence
window = 150 - 25 = 125 bp -> substring position 0-125, 20-145, 40 - 165, 60 - 185 of VR sequence
window = 125 - 25 = 100 bp -> substring position 0-100, 20-120, 40-140, 60-160, 80-180, 100-200
window = 100 - 25 =  75 bp -> search stop here coz window length (75bp) is less than min match length threshold -l 100bp

";

die $usage unless defined $opt_i and defined $opt_b;
die $usage . "\n" . $usage_long if defined $opt_H;
die "fastq file -i $LCY$fqFile$N doesn't exist!\n" if not -e $fqFile;
die "query file -b $LCY$VRFile$N doesn't exist!\n" if not -e $VRFile;
my $step_size            = defined $opt_s ? $opt_s : 25;
my $window_subtract_size = defined $opt_w ? $opt_w : 25;
my $minlen_threshold     = defined $opt_l ? $opt_l : 100;
my $outDir               = defined $opt_o ? $opt_o : "./";
my $quiet = 1 if defined $opt_q;
#
#" unless @ARGV;
my %VR;
system("mkdir -p $outDir") if not -d $outDir;
$outDir =~ s/\/+$//;
my ($folder1, $fileName1) = getFilename($fqFile, "folderfull");
my ($folder2, $VRFilename) = getFilename($VRFile, "folderfull");
my $logFile = "$outDir/$fileName1\_$VRFilename\_l$minlen_threshold\_w$window_subtract_size\_s$step_size.log";
my $outFile = "$outDir/$fileName1\_$VRFilename\_l$minlen_threshold\_w$window_subtract_size\_s$step_size.tsv";
my $outFileshort = "$outDir/$fileName1\_$VRFilename\_l$minlen_threshold\_w$window_subtract_size\_s$step_size.short.tsv";
my $VRstatsFile = "$outDir/$fileName1\_$VRFilename\_l$minlen_threshold\_w$window_subtract_size\_s$step_size.stats";

open (my $outLog, ">", $logFile) or die "Can't write to $logFile: $!\n";
my %out;
LOG($outLog, "$LGN
Parameters:$N
-i fastq file            : $LCY$fqFile$N
-b VR/query fasta file   : $LCY$VRFile$N
-o output dir            : $LCY$outDir$N
-s step size             : $LGN$step_size$N
-w window subtract size : $LGN$window_subtract_size$N
-l min length threshold  : $LGN$minlen_threshold$N

",$opt_q);
LOG($outLog, "\n#queries:\n",$opt_q);
open (my $in0, "<", "$VRFile") or die;#checker2_queries.tsv");
while (my $line = <$in0>) {
	chomp($line);
	my ($VRname, $VR);
	if ($line =~ /^>/) {
		($VRname) = $line =~ /^>(.+)$/;
		$line = <$in0>; chomp($line);
		($VR) = $line =~ /^(.+)$/;
	}
	else {
		($VRname, $VR) = split("\t", $line);
	}
	my $VR2 = $VR;
	$VR2 =~ s/^(.)/N/;
	$VR{$VR} = $VRname;
	#$VR{$VR2} = $VRname;
	my $VRshort = $VR;
	($VRshort) = $VR =~ /^(.{20})/ if length($VR) >= 20;
	print "$LCY$VRname$N\t$YW$VRshort$N...\n" if not defined $opt_q;
	LOG($outLog, "$VRname\t$VR\n",1);
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
	my @bcseq;
	#print "foreach\n" if ($read_name eq "\@6acb2984-5480-415b-96db-da1208fa1a0d");
	foreach my $VRorig (sort keys %VR) {
		my $minlen = 200;
		my $found = 0;
		my $VRrevorig = revcomp($VRorig);
		#print "   while $LGN$VR{$VRorig}$N\n" if ($read_name eq "\@6acb2984-5480-415b-96db-da1208fa1a0d");
		while ($minlen >= $minlen_threshold) {
		my $i = -1;
		while (1) {
			my ($VR, $VRrev);
			if ($i == -1) {
				($VR, $VRrev) = ($VRorig, $VRrevorig);
			}
			else {
				my $substr = $minlen > length($VRorig) ? length($VRorig) : $minlen;
				$VR = substr($VRorig, $i, $substr);
				$VRrev = revcomp($VR);
			}
			if ($read_seq =~ /($VR|$VRrev)/) {
				if ($read_seq =~ /$VR/) {
					($VRseqbeg, $VRseqend) = $read_seq =~ /^(.*)$VR(.*)$/;
					$VRfinal = $VR{$VRorig};
					$VRseqfinal = $VR;
					$VRstrand = "+";
				}
				elsif ($read_seq =~ /$VRrev/) {
					($VRseqbeg, $VRseqend) = $read_seq =~ /^(.*)$VRrev(.*)$/;
					$VRfinal = $VR{$VRorig};
					$VRseqfinal = $VRrev;
					$VRstrand = "-";
				}
				$totalbc ++;
				$beg = length($VRseqbeg);
				$end = length($read_seq) - length($VRseqend);
				my $VRseqbegdot = join("", ("-") x length($VRseqbeg));
				my $VRseqenddot = join("", ("-") x length($VRseqend));
				#if ($totalbc == 1) {
				push(@bc, $VRfinal);
				push(@bcseq, $VRseqfinal);
				#	$data{sample}{$VRfinal} ++;
				#	$data{samplebc}{$VRseqfinal} ++;
				#}
				if ($totalbc == 2) {
					$data{multiplebc} ++;
				}
				#my $length = $end - $beg;
				my $length = length($VRseqfinal);
				print $out0 "$read_name\t$beg\t$end\t$VRfinal\t$totalbc;$i;$length\t$VRstrand\t$VRseqbegdot$VRseqfinal$VRseqenddot\t$read_seq\t$read_qual\n";
				#print "$YW$read_name$N: $LCY$VRfinal$N ($beg-$end $VRstrand)\n$LGN$read_seq$N\n$LPR$VRseqbegdot$LGN$VRseqfinal$N$LPR$VRseqenddot$N\n\n";
				#print "   i=$LGN$i$N: LCY$VRfinal$N $LGN$VR{$VRorig}$N\n" if ($read_name eq "\@6acb2984-5480-415b-96db-da1208fa1a0d");
				$found = 1;
				last if $found eq 1;
			}
			$i += $step_size;
			last if length($VRorig) < $i + $minlen;
			last if $found eq 1;
			}
			$minlen -= $window_subtract_size;
			last if length($VRorig) < $i + $minlen;
			last if $found eq 1;
		}
		#print "  while done\n" if ($read_name eq "\@6acb2984-5480-415b-96db-da1208fa1a0d");
	}
	if ($VRfinal eq "UNKNOWN") {
		$data{sample}{none} ++;
		$VRseqfinal = join("", ("-") x length($read_seq));
		print $out0 "$read_name\t$beg\t$end\t$VRfinal\t0\t$VRstrand\t$VRseqfinal\t$read_seq\t$read_qual\n";
		print $out0b "$read_name\tUNKNOWN\n";
		#print "$YW$read_name$N: $LCY$VRfinal$N\n$LGN$read_seq$N\n";
	}
	else {
		$VRfinal = join(";", @bc);
		$VRseqfinal = join(";", @bcseq);
		$data{sample}{$VRfinal} ++;
		$data{samplebc}{$VRseqfinal} ++;
		print $out0b "$read_name\t" . join(";", @bc) . "\n";
	}
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


