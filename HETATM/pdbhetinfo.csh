#!/bin/csh -f

if($#argv != 1  ) then 
  echo "Usage : <list>"
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 
set list = $1
set EXEC = hetinfoextract.csh
setupCommandsFromListAndSchedule.pl -lis  $list -sleep 10 -inter 1000 -string "pdbhetinfoSingleInfo.pl -pr REPLACE" -out $EXEC

echo "pdbhetinfoCollate.pl  -lis $list" >> Sch.$EXEC
runInBack.csh Sch.$EXEC 

