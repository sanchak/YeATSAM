#!/bin/csh -f

source $SRC/GNM/step.setup.csh 

setenv FASTADIR FASTADIR_NT
# split these into individual fasta files
if( ! -e $LIST) then 
    uniquifyFasta.pl -in $TRSFILE -out trinity -write 1
endif 

## optimized for sleep anf interval
setupCommandsFromListAndSchedule.pl -lis  $LIST -sleep 10 -inter 1000 -string "blastn -db $WGS -query $FASTADIR/REPLACE.ALL.1.fasta -out BLASTOUT_2WGS/REPLACE.blast.nt " -out exec.wgs

setupCommandsFromListAndSchedule.pl -lis  $LIST -sleep 3 -inter 2000 -string "blastn -db  $RRNADB -query $FASTADIR/REPLACE.ALL.1.fasta -out BLASTOUT_RRNA/REPLACE.blast.nt" -out exec.rrna

orfprocessList.csh $LIST 

setupCommandsFromListAndSchedule.pl -lis list.transcriptome.clean.ORFS -out fff -sleep 3 -inter 5 -string \ 'blastp -db DB/plantpep.fasta -query FASTADIR_ORF/REPLACE.ALL.1.fasta -out BLASTOUT_AA/REPLACE.blast.nt '

parseBlastLatestList.pl -out RRNA.anno -lis list.transcriptome.clean -blastdir BLASTOUT_RRNA/ -blastcutoff 1000 -forWGS 1 -findcha 0  > & ! /dev/null & 

cat RRNA.anno.1000.anno.real list.transcriptome.contam > list.transcriptome.contam.andRRNA


processSplitTRSErrors.pl -in LISTS/finalanno.75.anno.real -outf LISTS/list.finalanno -dirin FASTADIR_ORF -diro FASTADIR_ORFBEST

### this is for getting the scaff2trs.500 - the highscoring mapping of trs to the WGS
parseBlastLatestList.pl -out LISTS/map.full.TRS2Scaffold -lis list.transcriptome.clean -blastdir BLASTOUT_2WGS/ -blastcutoff 500 -forWGS 1 > & ! /dev/null &
## still reversed - only about 70 
checkBLASTtrs2scaffolddirection.pl -blastdi BLASTOUT_SCAFFOLD -out blast.err 
## lot of hardcoded stuff here
parseBlastLatestpairwiseAllScaffolds.pl -blas BLASTOUT_SCAFFOLD -just 0 -inf scaff2trs.500

cat INFO/*.info > ! III  &
echo Run this genReportsForScaff.pl -outf uuuu -inf III -ann genome.ann.realgood -ignore ~/ooo 

mkdir  ORFTMP 
mkdir  ORF 
setupCommandsFromListAndSchedule.pl -lis $LIST -out fff -sleep 1 -inter 100 -string 'getorf FASTADIR_NT/REPLACE.ALL.1.fasta ORFTMP/REPLACE.orf' 

## this renames the strings, adding ORF
setupCommandsFromListAndSchedule.pl -lis $LIST -out fff -sleep 1 -inter 100 -string 'fixORFnames.pl -inf ORFTMP/REPLACE.orf -out ORF/REPLACE.orf '

## creates the top 3 ORFs in FASTADIR_ORF, and the longest in FASTADIR_JUSTONLONGEST
setupCommandsFromListAndSchedule.pl -lis $LIST -out fff -sleep 1 -inter 100 -string 'findrepeatedorfs.pl -trs REPLACE -orfdir ORF/ -write 4 '

$SRC/GNM/catANNintOnefile.csh FASTADIR_ORF list list.transcriptome.clean.ORFS # you need to edit this file


setupCommandsFromListAndSchedule.pl -lis list.transcriptome.clean.ORFS -out fff -sleep 3 -inter 5 -string \
 'blastp -db DB/plantpep.fasta -query FASTADIR_ORF/REPLACE.ALL.1.fasta -out BLASTOUT_AA/REPLACE.blast.nt '


setupCommandsFromListAndSchedule.pl -lis list.transcriptome.clean.ORFS -out fff -sleep 3 -inter 100 -string 'parseWebBlastOuts.pl -outf jjj -tr REPLACE -cutoff 1E-8 -blast BLASTOUT_AA -orfdir FASTADIR_ORF -ver 0 '

# Check for ANN/ERR and ANN/WARN (repeats) 
# then run
catANNintOnefile.csh ANN/COMMANDS c commandstocopy.csh
source commandstocopy.csh > & ! /dev/null &  # this creates all the FASTA files in FASTADIR_ORFBEST 

## Make each fasta file have the right name
mv FASTADIR_ORFBEST FASTADIR_ORFBESTRAW
setupCommandsFromListAndSchedule.pl -lis list.transcriptome.clean.expanded -out fff -sleep 2 -inter 100 -string 'fixFastaNames.pl -in FASTADIR_ORFBESTRAW/REPLACE.ALL.1.fasta -trs REPLACE -out FASTADIR_ORFBEST/REPLACE.ALL.1.fasta '





