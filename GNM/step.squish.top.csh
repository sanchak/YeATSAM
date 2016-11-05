set INITLIST = trinity.list.clean.BLAST.anno.0.merged.uniq.mfix.substring

step.squish.single.csh kmerExpandRetroPlant kmerMerged2Retro 20 10  DB/cores-database $INITLIST
step.squish.single.csh kmerExpandLRR kmerMerged2LRR 20 10 DB/DB.LRR.404.fasta $INITLIST




TW LISTS/kmerMerged2LRR/kmerMerged2LRR.75.anno.real.0 LISTS/kmerMerged2Retro/kmerMerged2Retro.75.anno.real.0 

cat LISTS/kmerMerged2LRR/kmerMerged2LRR.75.anno.real.0 LISTS/kmerMerged2Retro/kmerMerged2Retro.75.anno.real.0 > ! CCC
UNIQ CCC 
TW $INITLIST CCC.uniq
mv -f ofhinAbutnotinB $INITLIST.nonretroLRR


set SELF=$INITLIST.nonretroLRR
\rm -f DB/DB.$SELF
createDBFasta.csh $SELF FASTADIR_ORFBEST DB/DB.$SELF junkcomm


#set ksize = 20 
#testkmer.pl -genom  DB/DB.$SELF -out kmerMergedNonRetroLRR.Self  -lis $SELF -ksize $ksize -isNT 0 -fasta FASTADIR_ORFBEST -hopping 0 
#mappingAddCount.pl -in kmerMergedNonRetroLRR.Self/kmerMergedNonRetroLRR.Self.0.$ksize.mapall 
#groupBasedonCutoff.pl -in kmerMergedNonRetroLRR.Self/kmerMergedNonRetroLRR.Self.0.$ksize.mapall.mapped.pw -out kmerMergedNonRetroLRR.Self/group.$ksize -cutoff 0 
#cat kmerMergedNonRetroLRR.Self/group.$ksize.first kmerMergedNonRetroLRR.Self/kmerMergedNonRetroLRR.Self.0.$ksize.notmapped > ! kmerMergedNonRetroLRR.Self/geneset.$ksize
#
#wc -l kmerMergedNonRetroLRR.Self/geneset.$ksize
