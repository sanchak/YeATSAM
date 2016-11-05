#!/bin/csh  -f

if($#argv != 3   ) then 
    echo "Usage: <DBFILE> <cutoff>  <OUTFILE>"
	exit 
endif 

## Remember you can use this in parallel for DBFILE
set DBFILE = $1 
set cutoff = $2 
set OUTFILE = $3 

set BLASTDIR = $DBFILE.BO

mkdir -p $BLASTDIR


setenv FASTADIR $DBFILE.uniq.FDIR

if( ! -e $DBFILE.uniq.psi) then
    fixFastaNm.pl -in $DBFILE -out $DBFILE.tmp -same 
    uniquifyFasta.pl -in $DBFILE -out $DBFILE.uniq
	makeblastdbP.csh $DBFILE.uniq
endif

#if(! -e $DBFILE.uniq.FDIR) then 
splitFasta.pl -in $DBFILE.uniq
#endif 

set list = $DBFILE.uniq.list
set maplenfile = $DBFILE.uniq.length.sort


setupCommandsFromListAndSchedule.pl -lis $list -out blast.$OUTFILE.csh 
-sleep 10 -inter 100 -string "blastp -db $DBFILE.uniq -query 
$FASTADIR/REPLACE.ALL.1.fasta -out $BLASTDIR/REPLACE.blast.nt " 
 
echo "parseBlastLatestList.pl -out $OUTFILE -list $list -blastdir $BLASTDIR/ -blastcutof $cutoff -mapf $maplenfile -forWGS 0 -findcha 0 -isNT=0 "

