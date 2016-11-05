
# carried over
#createDBFasta.csh list.merged.overlap FASTADIR_ORFBEST DB/DBFastabestMergedOverlap

makeblastdbP.csh DB/DBFastabestMergedOverlap
mkdir BLASTOUT_SUBSET
setupCommandsFromListAndSchedule.pl -lis list.merged.overlap  -out fff -sleep 3 -inter 100 -string 'blastp -db DB/DBFastabestMergedOverlap -query FASTADIR_ORFBEST/$i.ALL.1.fasta -out BLASTOUT_SUBSET/$i.blast.nt '


# this creates the group.nr.50 
parseBlastLatestList.pl -out group.nr -lis list.merged.overlap -blastdir BLASTOUT_SUBSET/ -blastcutoff 50 -forW 0 
   rseBlastLatestList.pl -out group.nr -lis list.merged.overlap -blastdir BLASTOUT_SUBSET/ -blastcutoff 50 -forW 0
## ensure each pair has only one entry - the highest...
groupUniqify.pl -out group.nr.50.fix -inf group.nr.50

processGroups.pl -outf gene.grouphead -inf group.nr.50.fix -lis list.merged.overlap -map2len map.name2length -map2scaffold map.TRS2scaffold -ignorefile ign -final 125

