#!/bin/csh -f

if($#argv != 6  ) then 
  echo "<name> <list>  <BLASTOUT> <findchar> <isNT> <strict>"
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

set PWD = ` pwd`
set name = $1
set list = $2
set BLASTOUT = $3
set findchar = $4
set isNT = $5
set strict = $6

mkdir -p $name 
if(! -e $name/$name.50) then 
   parseBlastLatestList.pl -out $name/$name -lis $list -blastdir $BLASTOUT -blastcutoff 50 -forWGS 0 -findchar $findchar -isNT $isNT -str $strict
endif 


foreach blastcutoff ( 100 150 300)
      annConvertToBlastCutoff.pl -out POTR2MALUS/POTR2MALUS -blastcutoff $blastcutoff -oldb 50
end
