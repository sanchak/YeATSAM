#!/bin/csh -f
if($#argv != 1  ) then 
  echo "<pdb>"
  exit 
endif 

set PDB = $1 
set list=$PDBDIR/$PDB.list
#findMaxDistInDiffPDBPairwise.pl -list $list -out PDBINFO/$PDB.contactchains -cu 3.7 -ah PDBINFO/$PDB.AHBS -pairwise 1 -tag $PDB 
findMaxDistInDiffPDBPairwise.pl -list $list -out PDBINFO/$PDB.contactchains -cu 12 -ah PDBINFO/$PDB.AHBS -pairwise 1 -tag $PDB 

