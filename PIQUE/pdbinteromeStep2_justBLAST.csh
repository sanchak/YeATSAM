#!/bin/csh -f


mkdir -p WORKDIR


set DBFILE=DB/DBAllPDBSANDCHAINS
unlink $DBFILE
createDBFasta.csh list.final.uniq.wrote $FASTADIR/  $DBFILE

makeblastdbP.csh DB/DBAllPDBSANDCHAINS

newfile.csh runblast.csh 
foreach i  (` cat  list.final.uniq.wrote `)
	if(! -e BLASTOUT/$i.blast.nt) then 
         echo BLASTP $DBFILE $FASTADIR/$i.ALL.1.fasta BLASTOUT/$i.blast.nt >> runblast.csh 
	endif 
end
wc -l runblast.csh 

#scheduleprocessInsertingsleep.pl -inf runblast.csh -int 10 -sleep 2 
#runInBack.csh Sch.runblast.csh

