#!/bin/csh -f

##############################################################3
## This step takes the longest - running BLAST on 3*ORFS ...
## We also create a BLASTDB of all ORFs, just for searching later on
##
## Next, we will process these blast results
##############################################################3
source $SRC/GNM/step.setup.csh 

setenv FASTADIR FASTADIR_ORF
setupCommandsFromListAndSchedule.pl -lis  $ORFS3  -out fff -sleep 3 -inter 5 -string  "blastp -db $DB_REFPROTEOME -query $FASTADIR/REPLACE.ALL.1.fasta -out $BO_AA/REPLACE.blast.nt "



