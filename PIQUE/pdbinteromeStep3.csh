#!/bin/csh -f

if($#argv != 1  ) then 
  echo "Usage : <PERCENTHOMOLOGY>"
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "


mkdir -p WORKDIR

set listquery=list.final.uniq.wrote
set PERCENTHOMOLOGY=$1

set ignorepdbs=`cat ignorepdbs`


if(! -e mapfile.length) then 
   newfile.csh mapfile.length 
   foreach i  (` cat $listquery `)
       findfastalength.pl -inf $FASTADIR/$i.ALL.1.fasta >> mapfile.length
   end 
endif 

parseBlastLatestList.pl -out group.nr -lis $listquery -blastdir BLASTOUT/ -blastcutoff $PERCENTHOMOLOGY -forW 0 -findc 0 -stri 0 -isNT 0 -map mapfile.length

set GRPFILE = LISTS/group.nr/group.nr.$PERCENTHOMOLOGY

cat $GRPFILE | egrep -v "$ignorepdbs" >  ! JJJ
\mv -f JJJ WORKDIR/group.nr.$PERCENTHOMOLOGY

groupBasedonCutoff.pl -outf WORKDIR/finalgroup -inf $GRPFILE -cutoff $PERCENTHOMOLOGY -dir 1
grep only $GRPFILE > ! onlyone
extractindexfromfile.pl -in onlyone
cat onlyone.0 >> WORKDIR/finalgroup.first

## Now you need to finalgroup.first.ANNO by HAND 

\mv -f namemap namemap.old

cat onlyone | egrep -v "4X23|5BS7" >> finalgroup
UNIQ finalgroup 
#./pdbinteromeNameMap.pl -outf namemap -inf finalgroup -map finalgroup.first.ANNO
pdbinteromeNameMap.pl -outf WORKDIR/namemap -inf WORKDIR/finalgroup -map finalgroup.first.ANNO


