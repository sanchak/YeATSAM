#!/bin/csh -f


mkdir -p WORKDIR

set ignorepdbs=`cat ignorepdbs`


echo "We will ignore $ignorepdbs"

## HACK - remove this PDB its a chimeric one, fused
cat listofPDB5 | egrep -v "$ignorepdbs" > ! list.final
setenv FASTADIR FASTADIR
UNIQ list.final
writeFasta.csh list.final.uniq

set DBFILE=DB/DBAllPDBSANDCHAINS
unlink $DBFILE
createDBFasta.csh list.final.uniq.wrote $FASTADIR/  $DBFILE makeblastdbP.csh

newfile.csh runblast.csh 
foreach i  (` cat  list.final.uniq.wrote `)
	if(! -e BLASTOUT/$i.blast.nt) then 
         echo BLASTP $DBFILE $FASTADIR/$i.ALL.1.fasta BLASTOUT/$i.blast.nt >> runblast.csh 
	endif 
end
wc -l runblast.csh 

#scheduleprocessInsertingsleep.pl -inf runblast.csh -int 10 -sleep 2 
#runInBack.csh Sch.runblast.csh

