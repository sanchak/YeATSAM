#!/bin/csh -f

if($#argv != 3  ) then 
  echo "<TRS> <DB> <BLASTDIR>"
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

set PWD = ` pwd`
set TRS = $1
set DB = $2
set blastdir = $3


if(! -e $blastdir/$TRS.blast.nt  || -z $blastdir/$TRS.blast.nt) then 
    orfprocessSingle.csh $TRS
    BLASTP $DB FASTADIR_JUSTONLONGEST/$TRS.ALL.1.fasta $blastdir/$TRS.blast.nt 
    \rm -f  FASTADIR_JUSTONLONGEST/$TRS.ALL.1.fasta 
endif 
