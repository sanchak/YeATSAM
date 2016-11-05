#!/bin/csh -f

if($#argv != 3  ) then 
  echo "Usage :  <command> <inputfile> <output> : command can be zcat or  cat depending on whether files .gz"
  exit 
endif 


$1 $2 | awk 'NR%4==1{printf ">%s\n", substr($0,2)}NR%4==2{print}' > $3
