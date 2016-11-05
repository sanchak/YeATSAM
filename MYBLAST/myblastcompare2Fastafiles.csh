#!/bin/csh -f
if($#argv != 4  ) then 
  echo "Usage : <fasta1> <fasta2> <out> <what - P or N >"
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 

#mkdir -p ~/junk/TMPFASTADIR
#\cp -f $1 $2 ~/junk/TMPFASTADIR

 
myblastcompare2.pl -p1 $1 -p2 $2 -filesgiven -outf $3 -what $4

