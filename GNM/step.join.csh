#!/bin/csh -f

if($#argv != 3  ) then 
  echo "usage: <mergedir> <annfile> <DBfile>"
  exit 
endif 

set mergedir = $1
#set annfile = trinity.list.clean.BLAST.anno
set annfile = $2
set DBFILE = $3

mkdir -p $mergedir

if(! -e  $mergedir/joinscript.1) then 
   uniquifyFasta.pl -in $DBFILE
   joinTRS.csh $DBFILE.uniq  $mergedir
   source $mergedir/joinscript.1
   exit ;
endif 

cat $mergedir/*results > ! $mergedir/MERGEALL 
$SRC/GNM/joinTRSMergeList.pl -DB DB/plantpep.fasta -fas FASTADIR_ORFBEST -merge $mergedir -blast BLASTOUT_JOIN -out runjoin.csh -annfile $annfile


\rm -f $mergedir/newDB

echo Now run: createDBFasta.csh $mergedir/newlist.1  FASTADIR_ORFBEST/ $mergedir/newDB junkcomm


