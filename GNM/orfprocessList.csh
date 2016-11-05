#!/bin/csh -f

if($#argv != 3  ) then 
  echo "Usage <list>  <output>  <FDIR>"
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

set list = $1
set output = $2
set FDIR = $3

mkdir -p ORF
mkdir -p ORFTMP

echo "note we check FASTADIR_JUSTONLONGEST"

newfile.csh $output
foreach i (`cat $list`)
	if(! -e FASTADIR_JUSTONLONGEST/$i.ALL.1.fasta) then 
        echo orfprocessSingle.csh $i $FDIR >>  $output
	endif 
end 

scheduleprocessInsertingsleep.pl -int 300 -sleep 5 -inf $output
wc -l $output
#runInBack.csh Sch.ORFCommands.csh 
