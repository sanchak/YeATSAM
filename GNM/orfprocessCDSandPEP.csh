#!/bin/csh -f

if($#argv != 2  ) then 
  echo "Usage : <cds.fa>"
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

set TRSFILE=$1
setenv FASTADIR FASTADIR_NT
uniquifyFasta.pl -in $TRSFILE -out list.cds -write 1

set AAFILE=$2
setenv FASTADIR FASTADIR_AA
uniquifyFasta.pl -in $TRSFILE -out list.aa -write 1
TW list.cds.list list.aa.list 


orfOnlyList.csh list.cds.list
scheduleprocessInsertingsleep.pl -int 500 -sleep 20 -inf runallorf.csh
runInBack.csh Sch.runallorf.csh




