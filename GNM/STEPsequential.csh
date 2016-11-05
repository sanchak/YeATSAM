#!/bin/csh -
if($#argv != 3  ) then
  echo "Usage : <TAG> <trsfile> <DB2MATCH>"
  exit 
endif

set TAG=$1
set TRSFILE=$2
set DB2MATCH=$3

source $SRC/GNM/step.setup.csh $TAG
if(! -e DB) then 
    ln -s $WORKWALNUT/DB .
endif 

#\rm $LIST.* to force ...
setenv FASTADIR FASTADIR_NT
if( ! -e $LIST) then 
    uniquifyFasta.pl -in $TRSFILE -out $TAG -write 1
	extractFastaNames.pl -fast $TRSFILE -out $LIST
endif 

## optimized for sleep anf interval
if( -e  $DB_WGS) then 
   setupCommandsFromListAndSchedule.pl -lis  $LIST -sleep 10 -inter 1000 -string "blastn -db $DB_WGS -query $FASTADIR/REPLACE.ALL.1.fasta -out BLASTOUT_2WGS/REPLACE.blast.nt " -out exec.wgs
   set tag=trs2scaffold
   set blastdir=BLASTOUT_2WGS
   set blastcutoff=75
   parseBlastLatestList.pl -out $tag -lis $LIST -blastdir $blastdir -blastcutoff $blastcutoff -forWGS 1 -findcha 0 -isNT=1  
   \cp -f LISTS/$tag/$tag.$blastcutoff.anno.real.0  $LISTCLEAN
   source exec.wgs 
else
	echo " No WGS file $DB_WGS"
	\cp -f $LIST $LISTCLEAN
	sleep 1 
endif 


## ORFS 
if(! -e $ORFS3) then 
    getorf -sequence $TRSFILE -outseq $LIST.orf
    orfParseSingleFile.pl -out $ORFS3 -in $LIST.orf -cut 60
endif 



setenv FASTADIR FASTADIR_ORF
setupCommandsFromListAndSchedule.pl -lis  $ORFS3  -sleep 3 -inter 5 -string  "blastp -db $DB2MATCH -query $FASTADIR/REPLACE.ALL.1.fasta -out $BO_AA/REPLACE.blast.nt " -out exec.pep


#testkmer.pl -fastaf $TRSFILE -genom DB/viral.nt.fasta -out viral -ksiz 30 -isNT=1 -hopp 1 


set tag=allORFS3anno
set blastdir=$BO_AA
set blastcutoff=75
set list = $ORFS3
set expectedout = "LISTS/$tag/$tag.$blastcutoff.list"
echo parseBlastLatestList.pl -out $tag -lis $list -blastdir $blastdir -blastcutoff $blastcutoff -forWGS 1 -findcha 0 -isNT=1   >>  exec.pep 

echo "Running BLAST on ORFS..."
source exec.pep

## outfile is just dummy here
setupCommandsFromListAndSchedule.pl -lis  $LISTCLEAN -sleep 10 -inter 1000 -string "createANNFromBlastWithSplitORF.pl -tr REPLACE -blas BLASTOUT_AA -cutoff 1E-08   -outf ANN/COMMANDS/REPLACE.c" -out exec.anno -dontsleep 1
echo parseWebBlastOutsPostProcess.pl -lis $LISTCLEAN -tag $LISTCLEANBLAST >> exec.anno 
echo source $LISTCLEANBLAST.commands.csh >> exec.anno 
source  exec.anno 


