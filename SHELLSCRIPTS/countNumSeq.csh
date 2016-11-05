#!/bin/csh 

if($#argv != 1  ) then 
  echo "Usage : file"
  exit 
endif 

pruneSameSequenceFromMadeFasta.pl -outf tttttt -inf $1
 

