## this is for merging ORFS from one TRS with the same annotation
set INITLIST = trinity.list.clean.BLAST.anno.0.merged.uniq 

if( ! -e MERGEUNIQ/GOOD) then 
    scheduleprocessInsertingsleep.pl -inter 100 -sleep 5 -inf ANN/mergeUniqFixIssue.csh -out TMP/Sch.mergeuniq
    runInBack.csh TMP/Sch.mergeuniq
    cat MERGEUNIQ/*out | grep -v "BOTH" > ! MERGEUNIQ/GOOD
endif 

mergeUniqFixIssuePostProcess.pl -list $INITLIST -outf $INITLIST.mfix -inf MERGEUNIQ/GOOD
createDBFasta.csh $INITLIST.mfix FASTADIR_ORFBEST DB/DB.$INITLIST.mfix  makeblastdbP.csh 

SRC/GNM/getExactStrings.pl -outf $INITLIST.mfix.substring -fas DB/DB.$INITLIST.mfix

