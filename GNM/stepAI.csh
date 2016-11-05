#d!/bin/csh -f

source $SRC/GNM/step.setup.csh 


## create teh signal peptide DB from "known" MAKER annotations
if(! -e $DB_MAKERDB.SIGNALP) then 
    $SRC/GNM/genSignalPeptideFromDB.pl -outf $DB_MAKERDB.SIGNALP -in $DB_MAKERDB -size $SIZESIGNNALP
endif 

set AIINFO = $LISTCLEANBLASTNOANNO.info

## Certain proteins in $LISTCLEANBLASTNOANNO might not have start codon
## these are removed in LISTCLEANBLASTNOANNO_AI
if(! -e $AIINFO) then
   newfile.csh $AIINFO
   foreach i (` cat $LISTCLEANBLASTNOANNO` )
            cat $FA_AI/$i.info >> $AIINFO
   end
endif 

## Remove proteins < $CUTOFF_AI_SIZE
sort.pl -cutoff $CUTOFF_AI_SIZE -idx 1 -inf $AIINFO
extractindexfromfile.pl -in $AIINFO.sort

set LISTCLEANBLASTNOANNO_AI = $LISTCLEANBLASTNOANNO.info.sort.0



wc -l $LISTCLEANBLASTNOANNO*

if(! -e $DB_AI) then 
   createDBFasta.csh $LISTCLEANBLASTNOANNO_AI $FA_AI $DB_AI  makeblastdbP.csh
   $SRC/GNM/genSignalPeptideFromDB.pl -outf $DB_AI.SIGNALP -in $DB_AI -size $SIZESIGNNALP

   setenv FASTADIR FASTADIR_AI_SIGNALP
   splitFasta.pl -in $DB_AI.SIGNALP
endif 



setupCommandsFromListAndSchedule.pl -lis  $LISTCLEANBLASTNOANNO_AI -out exec.AI.csh -sleep 2 -inter 1000 -string  "blastp -db $DB_MAKERDB.SIGNALP -query FASTADIR_AI_SIGNALP/REPLACE.ALL.1.fasta -out $BO_AI/REPLACE.blast.nt " 

set blastcutoff = 35 
parseBlastLatestList.pl -out AI_SIGNALPMATCH -lis $LISTCLEANBLASTNOANNO_AI -blastdir $BO_AI -blastcutoff $blastcutoff -forWGS 1 -findcha 0 -isNT=1


## Uniquify - this should have been done earlier, but fine here too
createDBFasta.csh LISTS/AI_SIGNALPMATCH/AI_SIGNALPMATCH.$blastcutoff.anno.real.0 FASTADIR_ORFBEST.STARTCODON TMP/DDD junk
uniq
uniquifyFasta.pl -in TMP/DDD 

extractthoseinlist.pl -inf $AIINFO -list TMP/DDD.uniq.list -out AI 
aiNewProteins.pl -out AI.good -inf AI -cutoff $CUTOFF_AI_SIZE


sort.pl -in AI.good -idx 1
frequencyDistributionAbs.pl -outf AI.freq -inf AI.good.sort -max 2000 -delta 100 -start 0 -idx 1



## for paper tables
extractthoseinlist.pl -inf bwa_counts_run1.txt.0.CLEAN -list AI.good.0
countspreprocessBWA.pl -out jjj -in ttt.bwa_counts_run1.txt.0.CLEAN

sort.pl -in jjj -idx 2 ## sort on number of places
head -40 jjj.sort > ! H.jjj
sort.pl -in H.jjj -idx 1 -cut 400
extractthoseinlist.pl -inf bwa_counts_run1.txt.0.CLEAN -list H.jjj.sort.0 



