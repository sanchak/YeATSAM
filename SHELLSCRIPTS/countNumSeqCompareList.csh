#!/bin/csh 

if($#argv != 2  ) then 
  echo "Usage : file"
  exit 
endif 

pruneSameSequenceFromMadeFasta.pl -outf tttttt -inf $1
TW tttttt $2 
\rm tttttt
