#!/bin/csh -f

source $SRC/GNM/step.setup.csh 

## This is not of much use though - since we use createANNFromBlastWithSplitORF.pl to 
## analyze the outputs of the 3ORFS BLAST runs
set tag=allORFS3anno
set blastdir=$BO_AA
set blastcutoff=75
set list = $ORFS3
set expectedout = "LISTS/$tag/$tag.$blastcutoff.list"
if( ! -e $expectedout) then
   parseBlastLatestList.pl -out $tag -lis $list -blastdir $blastdir -blastcutoff $blastcutoff -forWGS 1 -findcha 0 -isNT=1  
endif 
## after running we might not get the output, so it should be a diff if
if( ! -e $expectedout) then
    echo "File expectedout: $expectedout not made, quitting"
	exit 
endif 


## outfile is just dummy here
setupCommandsFromListAndSchedule.pl -lis  $LISTCLEAN -sleep 10 -inter 1000 -string "createANNFromBlastWithSplitORF.pl -tr REPLACE -blas BLASTOUT_AA -cutoff 1E-08   -outf ANN/COMMANDS/REPLACE.c" -out exec.anno
echo parseWebBlastOutsPostProcess.pl -lis $LISTCLEAN -tag $LISTCLEANBLAST >> Sch.exec.anno 
echo scheduleprocessInsertingsleep.pl -inter 1000 -sleep 3 -inf $LISTCLEANBLAST.commands.csh >> Sch.exec.anno 
echo runInBack.csh Sch.$LISTCLEANBLAST.commands.csh >> Sch.exec.anno 
runInBack.csh Sch.exec.anno 
 

