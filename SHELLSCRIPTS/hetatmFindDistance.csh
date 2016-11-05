#!/bin/csh -f
if($#argv != 3  ) then 
  echo "Usage : <pdb> <het> <PDBDIR>"
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

mkdir -p TMPDIR 

set PDBDIR=$3

grep -v $2 $PDBDIR/$1.pdb > ! TMPDIR/$1.pdb
grep  $2 $PDBDIR/$1.pdb > ! TMPDIR/$2.pdb

setenv PDBDIR TMPDIR
findMaxDistInDiffPDB.pl  -p1 $1 -p2 $2
source $SRC/setup.csh 
