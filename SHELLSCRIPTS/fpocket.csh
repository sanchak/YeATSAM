#!/bin/csh -f

if($#argv != 1  ) then 
  echo "<PDB> "
  exit 
endif 


if(! -e $PDBDIR/$1_out/pockets/pocket0_atm.pdb) then 
   fpocket -f $PDBDIR/$1.pdb                                     
endif 

foreach l ( 0 1 2 3 4 )
     fpocketListResidues.pl -outf kkk -f $PDBDIR/$1_out/pockets/pocket${l}_atm.pdb -pr $1
end

