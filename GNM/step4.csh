#!/bin/csh -f

source $SRC/GNM/step.setup.csh 


## get info
setupCommandsFromListAndSchedule.pl -lis $LISTCLEANBLAST.list -sleep 2 -inter 100 -string "percentAminoAcid.pl  -inf FASTADIR_ORFBEST//REPLACE.ALL.1.fasta -pro REPLACE  -out FASTADIR_ORFBEST//REPLACE.info"  -out runstep3.csh 

## create FASTADIR_ORFBEST.STARTCODON for all - though we will need only for not annotated
setupCommandsFromListAndSchedule.pl -lis $LISTCLEANBLAST.list -sleep 2 -inter 100 -string "fastaSearchForStartCodon.pl -inf FASTADIR_ORFBEST/REPLACE.ALL.1.fasta -out FASTADIR_ORFBEST.STARTCODON/REPLACE.ALL.1.fasta" -out runstep3.csh -append
## and get infor for all too 
setupCommandsFromListAndSchedule.pl -lis $LISTCLEANBLAST.list -sleep 2 -inter 100 -string "percentAminoAcid.pl -inf FASTADIR_ORFBEST.STARTCODON//REPLACE.ALL.1.fasta -pro REPLACE -out FASTADIR_ORFBEST.STARTCODON//REPLACE.info" -out runstep3.csh -append




#processGroups.pl -outf gene.grouphead -inf group.nr.50.fix -lis list.merged.overlap -map2len map.name2length -map2scaffold map.TRS2scaffold -ignorefile ign -final 125

