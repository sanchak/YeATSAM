
createDBFasta.csh list.transcriptome.clean.expanded  FASTADIR_ORFBEST DB/DBFastabest

newfile.csh map.name2length 
pruneSameSequenceFromMadeFasta.pl -outf list.merged -inf DB/DBFastabest -length map.name2length


mappingAddCount.pl -out list.splices -inf list.merged.mapping -ignoresingle
createDBFasta.csh list.merged FASTADIR_ORFBEST DB/DBFastabestMerged


# this prints comands in joinscript.1  joinscript.2  
joinTRS.csh DB/DBFastabestMerged MERGEDALL100k

# these commands ensure that the merged TRS increase the balstscore when merged
setupCommandsFromListAndSchedule.pl -lis merge.list -out blastformerge -sleep 2 -inter 7 -string 'blastp -db DB/plantpep.fasta -query FASTADIR_ORFBEST/$i.ALL.1.fasta -out BLASTOUT_AA//$i.blast.nt '

setupCommandsFromListAndSchedule.pl -lis merge.list -out fff -sleep 30 -inter 1000 -string 'parsewebblast.pl -in BLASTOUT_AA//$i.blast.nt -outf BLASTOUT_AA/$i.ann'

joinTRSCheckBlast.pl -outf mergeresults -inf merge.list -blas BLASTOUT_AA/
TW list.merged merge.exclude
cat ofhinAbutnotinB merge.include > !  list.merged.overlap

joinTRSCreatePutativeList.pl -merge MERGEDALL100k -outf list.merge.putative


## Merging is a different step
#mergeEquivalent.pl -mapp list.splices -inf $countfile -out $countfile.merged

createDBFasta.csh list.merged.overlap FASTADIR_ORFBEST DB/DBFastabestMergedOverlap
