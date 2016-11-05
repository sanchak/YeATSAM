#!/bin/csh -f

if($#argv != 1  ) then 
  echo "Usage <PDBid>"
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

touch pymol.in
pymol.2oneprotein.pl -out $1.p1m -pdb1 $PDBDIR/$1.pdb -in pymol.in
pymol $1.p1m > & ! $1.pymol.log & 


