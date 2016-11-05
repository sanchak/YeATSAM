#!/bin/csh -f

if($#argv != 1  ) then 
  echo "Usage : <4 letter pdb id>"
  exit 
endif 


~/firefox/firefox-bin "http://www.rcsb.org/pdb/explore.do?structureId=$1"
