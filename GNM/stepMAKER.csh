#!/bin/csh -f

source $SRC/GNM/step.setup.csh 

set cutoff=60 


mkdir -p BLASTOUT_MAKER.$cutoff
setupCommandsFromListAndSchedule.pl -lis list.makermodels.list -sleep 2 -inter 50 -string "blastp -db $DB_ANNOuniqMERGED -query FASTADIR_MAKER/REPLACE.ALL.1.fasta -out BLASTOUT_MAKER.$cutoff/REPLACE.blast.nt" -out maker.csh 
setupCommandsFromListAndSchedule.pl -lis $LISTCLEANBLASTANNOMERGEDUNIQ -sleep 2 -inter 50 -string "blastp -db DB/walnut.wgs.5d.all.maker.proteins.fasta -query FASTADIR_ORFBEST/REPLACE.ALL.1.fasta -out BLASTOUT_MAKER.$cutoff/REPLACE.blast.nt " -out maker.csh  -append



