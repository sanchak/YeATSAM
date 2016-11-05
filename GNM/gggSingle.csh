#!/bin/csh -f
if($#argv != 1  ) then 
  echo "Usage : ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

set string2search = $1 
set dir = /home/sandeepc/DATA/GNM/


grep $string2search $dir/LISTS/annoGoodChar/annoGoodChar.60.anno 


echo " TRS   CE   CI   CK   EM   FL   HC   HL   HP   HU   IF   LE   LM   LY   PK   PL   PT   RT   SE   TZ   VB " 
grep $string2search $dir/COUNTS/bwa_counts_run1.txt.0.CLEAN.normalized  



grep $string2search $dir/ALLSCAFF.INFO

#wc -l grepres.*
