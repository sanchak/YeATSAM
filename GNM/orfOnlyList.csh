#!/bin/csh -f

if($#argv != 1  ) then 
  echo "Usage <list"
  exit 
endif 
mkdir -p ORF
mkdir -p ORFTMP
set list=$1 

newfile.csh runallorf.csh 
foreach i (`cat $list`)
   if(! -e ORF/$i.orf) then 
       echo getorf.csh $i >> runallorf.csh 
   endif 
end

