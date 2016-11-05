#!/bin/csh -f
if($#argv != 1  ) then 
  echo "<PDBid> "
  exit 
endif 


set i = $1 
### First move into PDBDIR - so set it properly
cd $PDBDIR 

set doit = 0 
if(! -e $i.pdb) then
   wget    http://www.rcsb.org/pdb/files/$i.pdb 
   set doit = 1
endif 


if( -e $PDBDIR/$i.list) then 
    foreach III (`cat $PDBDIR/$i.list`)
	    #echo "Searching for $PDBDIR/$III.pdb"
		if(! -e "$PDBDIR/$III.pdb") then 
			echo "$PDBDIR/$III.pdb does not exist"
			set doit = 1
		endif 
	end
else
	set doit = 1
endif 


if($doit == 1) then 
    splitpdbIntChains.pl -p $i  -outf $PDBDIR/$i.list
    #UNIQ $PDBDIR/$i.list
	getPDBModel1ChainAlist.csh $PDBDIR/$i.list
endif 

cd - 


