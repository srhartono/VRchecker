### Synopsis

Usage: ./VRchecker.pl -i <FqFile.fastq/fastq.gz> -b <VRFile.fa>

-o: output dir [currrent directory]

-l: min continuous match length threshold in bp [200]

-q: quiet, don't print anything

### Example:

- Default: Exact match `./VRchecker.pl -i example.fastq.gz -b example_VR.fa`

If run successfully it'll create files with the following md5sum:

```
372746bba70a0f1a18e77642dd0fc27b  example.fastq.gz_example_VR.fa_l0_w25_s20.log
cf5800bad17df36e63939ed37139e556  example.fastq.gz_example_VR.fa_l0_w25_s20.short.tsv
b73153444081faa831c08b5abe52e532  example.fastq.gz_example_VR.fa_l0_w25_s20.stats
13fe5ec4c35d783111140caf6c68945d  example.fastq.gz_example_VR.fa_l0_w25_s20.tsv
```

- At least 100bp of VR sequence match: `./VRchecker.pl -i example.fastq.gz -b example_VR.fa -l 100`

If run successfully it'll create files with the following md5sum:

```
7d469bbf7ece214c9ab81e660143aefc  example.fastq.gz_example_VR.fa_l100_w25_s20.log
3a3a19e3f10a79402f0911a926a1992b  example.fastq.gz_example_VR.fa_l100_w25_s20.short.tsv
95136045870286478649f3b2e4f6359f  example.fastq.gz_example_VR.fa_l100_w25_s20.stats
63765c04a0685d6c5aa4065023f3b93b  example.fastq.gz_example_VR.fa_l100_w25_s20.tsv
```
