#!/bin/csh -f

source $SRC/GNM/step.setup.csh 

set cutoff=60 

if(! -e  LISTS/wgs2maker) then 
   parseBlastLatestList.pl -out wgs2maker -list $LISTCLEANBLASTANNOuniqMERGED -blastdir BLASTOUT_MAKER.$cutoff/ -blastcutof $cutoff -forWGS 1 -findcha 0 -isNT=0
   parseBlastLatestList.pl -out maker2wgs -list list.makermodels.list -blastdir BLASTOUT_MAKER.$cutoff/ -blastcutof $cutoff -forWGS 1 -findcha 0 -isNT=0
endif 


set i=$cutoff 

\rm list.notin.*

set j=maker2wgs
\cp -f LISTS/$j/$j.$i.anno.no list.notin.YEATS
## remove retro
setupCommandsFromListAndSchedule.pl -lis list.notin.YEATS -sleep 10 -inter 500 -string "blastp -db  DB/DB.retroall -query FASTADIR_MAKER/REPLACE.ALL.1.fasta -out BLASTOUT_RETRO/REPLACE.blast.nt" -out UUUU
parseBlastLatestList.pl -out annoMakerRetro -list list.notin.YEATS -blastdir BLASTOUT_RETRO -blastcutoff 60 -forWGS 1 -findcha 1 -isNT=0
## find commin in plantsdb
parseBlastLatestList.pl -out maker2plants -lis LISTS/annoMakerRetro/annoMakerRetro.60.anno.no -blastdir BLASTOUT_MAKER.PLANTS/ -blastcutoff 60 -forWGS 1 -findcha 1 -isNT=0
\cp -f  LISTS/maker2plants/maker2plants.60.anno.real list.notin.YEATS.annot
excluderep.csh list.notin.YEATS.annot


set j=wgs2maker
TW LISTS/$j/$j.$i.anno LISTS/$j/$j.$i.anno.real
\mv  -f ofhinAbutnotinB list.notin.MAKER-P
setupCommandsFromListAndSchedule.pl -lis list.notin.MAKER-P -sleep 10 -inter 500 -string "blastp -db  DB/DB.retroall -query FASTADIR_ORFBEST/REPLACE.ALL.1.fasta -out BLASTOUT_RETRO/REPLACE.blast.nt" -out WWWW
 parseBlastLatestList.pl -out annoWalnutRetro -list list.notin.MAKER-P -blastdir BLASTOUT_RETRO -blastcutoff 60 -forWGS 1 -findcha 1 -isNT=0
extractthoseinlist.pl -inf LISTS/annoGoodChar/annoGoodChar.60.anno -list LISTS/annoWalnutRetro/annoWalnutRetro.60.anno.no -out list.notin.MAKER-P.noretro.annot
$SRC/GNM/excluderep.csh list.notin.MAKER-P.noretro.annot

exit 



####################
## Make NR set #####
##################

## takes idx one which includes the ORF
extractindexfromfile.pl -in list.notin.MAKER-P.annot -idx 1 
yeats_MakeNRSet.csh list.notin.MAKER-P.annot.1 DB/DB.list.notinMAKER  FASTADIR_ORF notin.MAKER-P 60 
parseBlastLatestList.pl -out BLASTOUTSUBSET/notin.MAKER-P.nonredundant -lis list.notin.MAKER-P.annot.1 -blastdir BLASTOUTSUBSET/ -blastcutof $cutoff -forWGS 0 -findcha 0 -isNT=1

extractindexfromfile.pl -in list.notin.YEATS -idx 1 
yeats_MakeNRSet.csh list.notin.YEATS.0 DB/DB.list.notinYEATS  FASTADIR_MAKER notin.YEATS 60 
parseBlastLatestList.pl -out BLASTOUTSUBSET/notin.YEATS.nonredundant -list list.notin.YEATS.0 -blastdir BLASTOUTSUBSET/ -blastcutof $cutoff -forWGS 0 -findcha 0 -isNT=1



##################################################
### Get Annotations for these grouped nonredundant 
##################################################


extractthoseinlist.pl -inf MAPDIRS/maker2PLANTS.60.anno.real -list BLASTOUTSUBSET/notin.YEATS.nonredundant.60.GROUPED -out BLASTOUTSUBSET/notin.YEATS.nonredundant.60.GROUPED.anno
 extractthoseinlist.pl -inf MAPDIRS/finalanno.60.anno.real -list BLASTOUTSUBSET/notin.MAKER-P.nonredundant.60.GROUPED -out BLASTOUTSUBSET/notin.MAKER-P.nonredundant.60.GROUPED.anno


