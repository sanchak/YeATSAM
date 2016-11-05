#!/bin/csh 

if($#argv != 3  ) then 
  echo "Usage : <pdb1>  <pdb2>  <outfile>"
  exit 
endif 


echo /home/sandeepc/DATA/externaltools/probis-2.4.2/probis -f1 $PDBDIR/$1.pdb -c1 A -f2 $PDBDIR/$2.pdb -c2 A -compare -super 
#/home/sandeepc/DATA/externaltools/probis-2.4.2/probis -f1 $PDBDIR/$1.pdb -c1 A -f2 $PDBDIR/$2.pdb -c2 A -compare -super 

parseProbis.pl -p1 $1 -p2 $2 -out $3

