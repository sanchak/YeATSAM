


If you have two fastq files your command:
bwa mem -M -t 16 ref.fa read1.fq read2.fq > aln.sam
is absolutely fine. Your read2.fq is called mates.fq in the bwa examples. If you view the first lines of both files the read names are identical despite a 1 or 2 for the corresponding read of the pair.

If you only have one interleaved fastq file you would use the -p option:
bwa mem -M -t 16 -p ref.fa read.fq > aln.sam
In this case both reads of a pair are in the same fastq file successively. Have a look at the read names.

For the unlikely case you would like to handle your paired-end reads as single ends the command is:
bwa mem -M -t 16 ref.fa read.fq > aln.sam



########## DO NOT PIPE STDERRR #############
bwa mem genome.fa  ED-N-1_CGGAAT_L002_R1_001.fastq.gz.fa ED-N-1_CGGAAT_L002_R2_001.fastq.gz.fa > ! ED-N-1.sam &


samtools view -f4 -S $i.sam > ! $i.unmapped.sam


bwa index $i.unmapped.fasta
