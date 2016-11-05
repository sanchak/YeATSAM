#!/bin/csh  -f

if($#argv != 2  ) then 
    echo "<inputfile> <list>"
endif 


extractthoseinlist.pl -inf $1 -lis $2
 

