#!/bin/csh -f
if($#argv != 2  ) then
  echo "Usage : <TAG> <trsfile>"
  exit 
endif

set TAG=$1
set TRSFILE=$2

source $SRC/GNM/step.setup.csh $TAG

setenv FASTADIR FASTADIR_NT
# split these into individual fasta files
if( ! -e $LIST) then 
    uniquifyFasta.pl -in $TRSFILE -out $TAG -write 1
	extractFastaNames.pl -fast $TRSFILE -out $LIST
endif 

## optimized for sleep anf interval
setupCommandsFromListAndSchedule.pl -lis  $LIST -sleep 10 -inter 1000 -string "blastn -db $DB_WGS -query $FASTADIR/REPLACE.ALL.1.fasta -out BLASTOUT_2WGS/REPLACE.blast.nt " -out exec.wgs

setupCommandsFromListAndSchedule.pl -lis $LIST -out exec.viral -sleep 10 -inter 1000 -string "blastn -db DB/viral.nt.fasta -query FASTADIR_NT/REPLACE.ALL.1.fasta -out BLASTOUT_VIRAL/REPLACE.blast.nt "

setupCommandsFromListAndSchedule.pl -lis  $LIST -sleep 3 -inter 2000 -string "blastn -db  $DB_RRNADB -query $FASTADIR/REPLACE.ALL.1.fasta -out BLASTOUT_RRNA/REPLACE.blast.nt" -out exec.rrna

## ORFS 
if(! -e $ORFS3) then 
    getorf -sequence $TRSFILE -outseq $LIST.orf
    orfParseSingleFile.pl -out $ORFS3 -in $LIST.orf -cut 60
endif 

# hack
\cp -f $TAG.list $TAG.list.clean

