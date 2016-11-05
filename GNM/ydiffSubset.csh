#!/bin/csh -f

if($#argv != 7  ) then 
  echo "Usage :  ydiffSubset <fastaAList> <fastaDirA>   <fastaBList>  <fastaBDIR> <comm> <cutoff> <tag>"
  echo "HEre A is a list now, with split in fastaDirA. fastaB is still a single file that will be made into a DB"
  echo "Makes sense to have lesser number in A"
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "


## TODO - cehck A has lesser number of fastas
set fastaALIST = $1
set fastaADIR = $2
set fastaBList = $3
set fastaBDIR = $4
set blast = $5
set cutoff = $6
set tag = $7

mkdir -p $tag 

set NEWDB = $tag/DB_B

createDBFasta.csh $fastaBList  $fastaBDIR  $NEWDB  # no need to make DB 



ydiffwithfDIR.csh $1 $2 $NEWDB $blast $cutoff  $tag
