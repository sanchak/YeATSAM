#!/bin/csh  -f

if($#argv != 3  ) then 
    echo "<filename> <percentchange> <overwrite> "
	exit 
endif 

set overwrite=$3

echo Writing to $2.$1
convert $1 -resize $2% $2.$1

if($overwrite == 1) then
  echo Overwriting $1 
  mv -f $2.$1 $1
endif 

