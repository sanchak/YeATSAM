#!/bin/csh -f
if($#argv != 1  ) then 
  echo "<pdb>"
  exit 
endif 

mkdir -p PDBINFO

set PDB = $1 

pdbgetlistAndSplit.csh $PDB

pdbAPBS.csh $PDB 

pdbgetPWDistance.csh $PDB
