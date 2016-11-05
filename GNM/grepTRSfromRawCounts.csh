#!/bin/csh  -f

if($#argv != 3  ) then 
  echo "<list> <rawcm> <out>"
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

set rawcount = $1
set list = $2 

newfile.csh $3
extractthoseinlist.pl -inf $rawcount -list $list  -tag ttt -inor
echo " TRS   CE   CI   CK   EM   FL   HC   HL   HP   HU   IF   LE   LM   LY   PK   PL   PT   RT   SE   TZ   VB " >> $3
cat ttt.$rawcount >> $3

 createTexTable.pl -in $3 
