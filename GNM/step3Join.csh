#!/bin/csh -f

source $SRC/GNM/step.setup.csh 

if( ! -e MERGEDIR.1/newDB) then 
   s $SRC/GNM/step.join.csh MERGEDIR $ANNFILE $DB_ANNOuniq
   s $SRC/GNM/step.join.csh MERGEDIR.1 MERGEDIR/newannfile.1 MERGEDIR/newDB
endif 


\cp -f MERGEDIR.1/newDB DB/DB.$LISTCLEANBLASTANNOMERGED
uniquifyFasta.pl -in DB/DB.$LISTCLEANBLASTANNOMERGED
make
\cp -f DB/DB.$LISTCLEANBLASTANNOMERGEDUNIQ.list $LISTCLEANBLASTANNOMERGEDUNIQ


## NOte findchar=1
parseBlastLatestList.pl -out annoGoodChar  -lis $ANNFILE.0  -blastdir BLASTOUT_AA/ -blastcutoff 60 -forWGS 1 -findcha 1 -isNT=0
ln -s LISTS/annoGoodChar/annoGoodChar.60.anno $ANNFILE.goodchar.anno


## Get merged annotations
parseBlastLatestList.pl -out annoMerge -list trinity.list.clean.BLAST.anno.mergeends.uniq -blastdir BLASTOUT_JOIN/ -blastcutoff 60 -forWGS 1 -findcha 1 -isNT=0 -strict 0
cat LISTS/annoMerge/annoMerge.60.anno.real $ANNFILE.goodchar.anno > ! $ANNFILE.goodchar.merge.anno
