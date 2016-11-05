#!/bin/csh  -f

if($#argv != 1  ) then 
    echo "Wrong args, required 1 " ; exit 
endif 

set list = $1 


source $list > & ! /dev/null &
 
set file = ~/lastrunprocess
touch $file 

echo #source $list >> $file
echo exit >> $file
echo kill -9 $! >> $file
reverseorder.pl -in $file
