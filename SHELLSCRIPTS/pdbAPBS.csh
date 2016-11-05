#!/bin/csh -f
if($#argv != 1  ) then 
  echo "<pdb>"
  exit 
endif 

set PDB = $1 

set list=$PDBDIR/$PDB.list

dssplist.csh  $list 
newfile.csh PDBINFO/$PDB.AHBS 
foreach i (`cat $list`)
   helixparseDSSPoutput.pl -outf PDBINFO/$PDB.AHBS -p $i -dssp $DSSP/$i.dssp -what HELIX -writeind 1
   #echo $status 
   helixparseDSSPoutput.pl -outf PDBINFO/$PDB.AHBS -p $i -dssp $DSSP/$i.dssp -what BETA -writeind 1
end

