#!/bin/csh -f

if($#argv != 1  ) then 
  echo "Usage : <cds.fa>"
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

set TRSFIL=$1
setenv FASTADIR FASTADIR_NT
uniquifyFasta.pl -in $TRSFILE -out list.cds -write 1
