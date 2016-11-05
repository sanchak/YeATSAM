#!/bin/csh -f

if($#argv != 6  ) then 
  echo "Usage :  ydiffbasicLIST <fastaAList> <fastaDirA>   <fastaB> <comm> <cutoff> <tag>"
  echo "HEre A is a list now, with split in fastaDirA. fastaB is still a single file that will be made into a DB"
  echo "Makes sense to have lesser number in A"
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "


## TODO - cehck A has lesser number of fastas
set fastaALIST = $1
set fastaADIR = $2
set fastaB = $3
set blast = $4
set cutoff = $5
set tag = $6


mkdir -p $tag 

ydiffwithfDIR.csh $1 $2 $3 $4 $5 $6
