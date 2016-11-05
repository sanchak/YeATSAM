#!/bin/csh  -f

if($#argv != 2  ) then 
    echo "Wrong args, required 1 " ; exit 
endif 

\mv -f "$1" ~/Documents/docsaved/$2 

 

