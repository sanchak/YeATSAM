#!/bin/csh -f

if($#argv != 1  ) then 
  echo "Usage : ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 

set i=$1 
\rm -rf velvet.$i
~/velvet_1.2.10/velveth velvet.$i 31 $i
 ~/velvet_1.2.10/velvetg velvet.$i  
