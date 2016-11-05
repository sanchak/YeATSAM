#!/bin/csh -f

if($#argv != 8  ) then 
  echo "Usage: <list>  <DB> <BLASTOUT> <EXEC> <FASTADIR> <FORCE> <INT> <SLEEP>"
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

set list = $1 
set DB = $2 
set BLASTOUT = $3 
set EXEC = $4 
set FASTADIR = $5
set FORCE = $6
set INT = $7
set SLEEP = $8


mkdir -p $BLASTOUT

newfile.csh $list.comm
foreach i (`cat $list`)
	 if(! -e $BLASTOUT/$i.blast.nt || $FORCE == 1) then 
         echo $EXEC $DB $FASTADIR/$i.ALL.1.fasta $BLASTOUT/$i.blast.nt  >> $list.comm
	 endif 
end 

wc -l $list
wc -l $list.comm


scheduleprocessInsertingsleep.pl -int $INT -sleep $SLEEP -inf $list.comm


echo "parseBlastLatestList.pl -out MMM -lis $list -blastdir $BLASTOUT -blastcutoff 300 -forWGS 0 -findcha 0 -isNT 1 -str 0"
