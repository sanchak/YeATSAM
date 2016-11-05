#!/bin/csh -f
if($#argv != 1  ) then 
  echo "Usage : ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

set list = $1 
set dir = /home/sandeepc/DATA/GNM/

extractthoseinlist.pl -lis $list -inf $dir/COUNTS/bwa_counts_run1.txt.0.CLEAN.normalized -ignorep -out $list.counts
echo " TRS   CE   CI   CK   EM   FL   HC   HL   HP   HU   IF   LE   LM   LY   PK   PL   PT   RT   SE   TZ   VB " >> $list.counts


extractthoseinlist.pl -lis $list -inf $dir/LISTS/annoGoodChar/annoGoodChar.60.anno -ignorep -out $list.anno
extractthoseinlist.pl -lis $list -inf $dir/ALLSCAFF.INFO -ignorep -out $list.scaffinfo


#grep $string2search $dir/LISTS/annoGoodChar/annoGoodChar.60.anno 
#
#grep $string2search $dir/COUNTS/bwa_counts_run1.txt.0.CLEAN.normalized  

#wc -l grepres.*
