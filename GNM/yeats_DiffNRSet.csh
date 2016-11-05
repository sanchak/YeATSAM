#!/bin/csh  -f

if($#argv != 7  ) then 
    echo "Usage: <list1> <FASTADIR1>  <list2> <FASTADIR2> <BLASTDIR> <FINALOUT> <cutoff>"
	exit 
endif 

## Remember you can use this in parallel for DBFILE
set list1 = $1 
set FASTADIR1 = $2 
set list2 = $3
set FASTADIR2 = $4 
set BLASTDIR = $5
set FINALOUT = $6 
set cutoff = $7 


set DBFILE2 = $list2.DBTMP
if(! -e $DBFILE2) then 
   createDBFasta.csh $list2 $FASTADIR2 $DBFILE2
   makeblastdbP.csh $DBFILE2
endif 


mkdir -p $BLASTDIR


echo ""
echo setupCommandsFromListAndSchedule.pl -lis $list1 -out blast.exec.csh -sleep 2 -inter 100 -string REPLACETICK blastp -db $DBFILE2 -query $FASTADIR1/\$i.ALL.1.fasta -out $BLASTDIR/\$i.blast.nt REPLACETICK 
echo ""
 
echo "parseBlastLatestList.pl -out $BLASTDIR///$FINALOUT -list $list1 -blastdir $BLASTDIR/ -blastcutof $cutoff -forWGS 0 -findcha 0 -isNT=0 "

