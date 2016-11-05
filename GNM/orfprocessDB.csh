#!/bin/csh -f

if($#argv != 2  ) then 
  echo "Usage <DB>  <output> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

set DBFILE = $1
set output = $2


setenv FASTADIR $DBFILE.FDIR

splitFasta.pl -in $DBFILE
extractFastaNames.pl -fas $DBFILE -out $DBFILE.list
orfprocessList.csh $DBFILE.list $2 $FASTADIR 
