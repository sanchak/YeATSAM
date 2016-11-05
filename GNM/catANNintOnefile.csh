#!/bin/csh 

if($#argv != 3  ) then 
  echo "<dir> <postfix> <outfile>"
  exit 
endif 
newfile.csh $3 
cat $1/C1*$2 >> $3
cat $1/C2*$2 >> $3
cat $1/C3*$2 >> $3
cat $1/C4*$2 >> $3

cat $1/C50*$2 >> $3
cat $1/C51*$2 >> $3
cat $1/C52*$2 >> $3
cat $1/C53*$2 >> $3
cat $1/C54*$2 >> $3
cat $1/C55*$2 >> $3
cat $1/C56*$2 >> $3
cat $1/C57*$2 >> $3
cat $1/C58*$2 >> $3
cat $1/C59*$2 >> $3

cat $1/C6*$2 >> $3
cat $1/C7*$2 >> $3
cat $1/C8*$2 >> $3
cat $1/C9*$2 >> $3
wc -l $3

