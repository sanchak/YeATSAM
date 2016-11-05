#!/bin/csh -f

if($#argv != 2  ) then 
  echo "<list of chr> <significance>"
  exit 
endif 

mkdir -p PARSEBLASTOUT 
set list=$1
set significance=$2
foreach chromosome (`cat $list`)
     parseBlastLatestList.pl -out PARSEBLASTOUT/MMM.$chromosome -lis LISTS/list.$chromosome -blastdir BLASTOUT -blastcutoff $significance -forWGS 0 -findcha 0 -isNT 1 -str 1
end


