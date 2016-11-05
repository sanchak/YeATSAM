### NEED TO FINS A WAY of gettign vars indisde the ticks
setupCommandsFromListAndSchedule.pl -lis list.transcriptome.clean.notinDB  -out startM.csh -sleep 2 -inter 100 -string ' fastaSearchForStartCodon.pl -inf FASTADIR_JUSTONLONGEST/$i.ALL.1.fasta -out FASTADIR_JUSTONLONGEST.STARTCODON/$i.ALL.1.fasta'

set MAKERDB=DB/walnut.wgs.5d.all.maker.proteins.fasta
set FASTADIRHERE=FASTADIR_JUSTONLONGEST.STARTCODON
set INITLIST=list.transcriptome.clean.notinDB.info
set AIDB = DB/DB.$INITLIST
set SIZE=20
set OUTPERCENTGOOD=$INITLIST.signalp.$SIZE.annovalues.sort.PERCENTGOOD


createDBFasta.csh $INITLIST.0 $FASTADIRHERE/ $AIDB

$SRC/GNM/genSignalPeptideFromDB.pl -outf $MAKERDB.SIGNALP -in $MAKERDB -size $SIZE
$SRC/GNM/genSignalPeptideFromDB.pl -outf $AIDB.SIGNALP -in $AIDB -size $SIZE
$SRC/GNM/compareSignalPeptideFromDB.pl -compa $MAKERDB.SIGNALP.$SIZE -in $AIDB -size $SIZE -out  $INITLIST.0.signalp.$SIZE

extractindexfromfile.pl -in $INITLIST.0.signalp.$SIZE
extractthoseinlist.pl -inf $INITLIST -list $INITLIST.0.signalp.$SIZE.0 -out $INITLIST.0.signalp.$SIZE.0.annovalues
sort.pl -in $INITLIST.0.signalp.$SIZE.0.annovalues -idx 1 -cutoff 100


$SRC/GNM/aiNewProteins.pl -inf $INITLIST.0.signalp.$SIZE.0.annovalues.sort -cutoff 100 -out $OUTPERCENTGOOD

createDBFasta.csh $OUTPERCENTGOOD $FASTADIRHERE/ $OUTPERCENTGOOD.DBTMP

makeblastdbP.csh $OUTPERCENTGOOD.DBTMP

yeats_MakeNRSet.csh $OUTPERCENTGOOD  $OUTPERCENTGOOD.DBTMP $FASTADIRHERE/ BLASTOUT_AI.$SIZE proteinAI.$SIZE 60

