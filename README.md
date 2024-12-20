### Synopsis

Usage: ./VRchecker.pl -i <File.fastq/fastq.gz> -b <VR.fa>

-o: output dir [currrent directory]
-l: min continuous match length threshold in bp [100]


### Example:

- Exact match only (assuming VR sequence length is 200bp)

`./VRchecker.pl -i 001_pFC8-VR-TAC-2x-LAC-Mini-1_reads.fastq.gz -b VR.fa -l 200`

- At least 100bp of VR sequence match:
  
`./VRchecker.pl -i 001_pFC8-VR-TAC-2x-LAC-Mini-1_reads.fastq.gz -b VR.fa -l 100`

