#!/bin/csh -f

if($#argv != 6  ) then 
  echo "<expanded name> <outdir> <cutoffhigh> <cutofflow> <DB2Expand>"
  exit 
endif 

set what = $1
set outdir = $2
set cutoffhigh = $3 
set cutofflow = $4 
set DB2Expand = $5
set list2squish = $6

## This expands the DB to include other proteins from DB/plantpep.fasta 

expandDBwithDB.csh $DB2Expand DB/plantpep.fasta FASTADIR_PLANTPEP/ $cutoffhigh $what

## This compares the merged set to the new expanded DB - at a lower significance, since we will
## make it higher with a pairwise.
## This slides on the genome too since the genome is small.
testkmer.pl -genom $what/DB.$what -out $outdir -lis $list2squish -ksize $cutofflow -isNT 0 -fasta FASTADIR_ORFBEST -hopping 0 


## However, this needs to be validated at a higher cutoff - so BLAST these
$SRC/MYBLAST/myblastcompare2FastafilesList.pl -out $outdir/comm -in $outdir/$outdir.0.$cutofflow.mapone -blast $what -f1 FASTADIR_ORFBEST -f2 FASTADIR_PLANTPEP/
scheduleprocessInsertingsleep.pl -inter 100 -sleep 3 -inf $outdir/comm -out $outdir/Sch.comm 
echo parseBlastLatestList.pl -out  $outdir -li $outdir/$outdir.0.$cutofflow.list -blastdir $what/ -blastcutoff 75 -forWGS 1 -findcha 0 -isNT=0  >> $outdir/Sch.comm 
runInBack.csh $outdir/Sch.comm



