#!/bin/csh -f

if($#argv != 1  ) then 
  echo "Usage : protein name "
  exit 
endif 


getPDBModel1ChainA.pl -protein $1 -model



#\rm $1*-m*.pdb $1*-c*.pdb 

#\rm $1.pdb
#mv -f A.pdb $1.pdb 
