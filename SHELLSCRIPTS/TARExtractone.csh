#!/bin/csh  -f

if($#argv != 2  ) then 
    echo "Wrong args, required 1 " ; exit 
endif 


echo tar -xf $1 --wildcards --no-anchored \'$2\*\'
 

