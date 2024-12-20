### Synopsis

Usage: ./VRchecker.pl -i <File.fastq/fastq.gz> -b <VR.fa>

-o: output dir [currrent directory]

-l: min continuous match length threshold in bp [100]


### Example:

- Exact match only (assuming VR sequence length is 200bp)

`./VRchecker.pl -i example.fastq.gz -b example_VR.fa -l 200`

- At least 100bp of VR sequence match:
  
`./VRchecker.pl -i example.fastq.gz -b example_VR.fa -l 100`

If run successfully it'll create files with the following md5sum:

- 12a01f26616ae2ef5242b8dcdb5ba572  example.fastq.gz_example_VR.fa_l200_w25_s20.log
- cf5800bad17df36e63939ed37139e556  example.fastq.gz_example_VR.fa_l200_w25_s20.short.tsv
- b73153444081faa831c08b5abe52e532  example.fastq.gz_example_VR.fa_l200_w25_s20.stats
- 13fe5ec4c35d783111140caf6c68945d  example.fastq.gz_example_VR.fa_l200_w25_s20.tsv
