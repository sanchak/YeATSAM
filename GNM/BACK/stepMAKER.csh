set cutoff=60 
\rm DB/DB.ORFBEST.$cutoff 
parseBlastLatestList.pl -out MAPDIRS/finalanno -lis list.transcriptome.clean.ORFS -blastdir BLASTOUT_AA/ -blastcutof $cutoff -forWGS 1 -findcha 0 -isNT=1
mkdir BLASTOUT_MAKER.$cutoff
processSplitTRSErrors.pl -in MAPDIRS/finalanno.$cutoff.anno.real -dirin FASTADIR_ORF -diro FASTADIR_ORFBEST
s MAPDIRS/finalanno.$cutoff.anno.real.parsed.commands.csh
createDBFasta.csh MAPDIRS/finalanno.$cutoff.anno.real.parsed.list FASTADIR_ORFBEST/ DB/DB.ORFBEST.$cutoff
makeblastdbP.csh DB/DB.ORFBEST.$cutoff


foreach i (`cat list.makermodels.list `)
  echo blastp -db DB/DB.ORFBEST.$cutoff -query FASTADIR_MAKER/$i.ALL.1.fasta -out BLASTOUT_MAKER.$cutoff/$i.blast.nt >> hhh
end

foreach i (`cat MAPDIRS/finalanno.$cutoff.anno.real.parsed.list `)
     echo blastp -db DB/walnut.wgs.5d.all.maker.proteins.fasta -query FASTADIR_ORFBEST/$i.ALL.1.fasta -out BLASTOUT_MAKER.$cutoff/$i.blast.nt  >> hhh
end 
scheduleprocessInsertingsleep.pl -int 500 -sleep 10 -inf hhh

parseBlastLatestList.pl -out MAPDIRS/wgs2maker -list MAPDIRS/finalanno.$cutoff.anno.real.parsed.list -blastdir BLASTOUT_MAKER.$cutoff/ -blastcutof $cutoff -forWGS 1 -findcha 0 -isNT=1
parseBlastLatestList.pl -out MAPDIRS/maker2wgs -list list.makermodels.list -blastdir BLASTOUT_MAKER.$cutoff/ -blastcutof $cutoff -forWGS 1 -findcha 0 -isNT=0


set j=maker2wgs
set i=60 
TW MAPDIRS/$j.$i.anno MAPDIRS/$j.$i.anno.real
\mv -f ofhinAbutnotinB list.notin.YEATS

set j=wgs2maker
set i=60 
TW MAPDIRS/$j.$i.anno MAPDIRS/$j.$i.anno.real
\mv  -f ofhinAbutnotinB list.notin.MAKER-P


## extractthoseinlist wont work for YEATS since we need ORF
mapANDextractthoseinlist.pl -map MAPDIRS/finalanno.60.anno.real.parsed.mapnames -inf MAPDIRS/finalanno.60.anno.real \
-lis MAPDIRS/wgs2maker.60.anno.no -out list.notin.MAKER-P.annot 
$SRC/GNM/excluderep.csh list.notin.MAKER-P.annot 


####################
## Make NR set #####
##################

## takes idx one which includes the ORF
extractindexfromfile.pl -in list.notin.MAKER-P.annot -idx 1 
yeats_MakeNRSet.csh list.notin.MAKER-P.annot.1 DB/DB.list.notinMAKER  FASTADIR_ORF notin.MAKER-P 60 
parseBlastLatestList.pl -out BLASTOUTSUBSET/notin.MAKER-P.nonredundant -lis list.notin.MAKER-P.annot.1 -blastdir BLASTOUTSUBSET/ -blastcutof $cutoff -forWGS 0 -findcha 0 -isNT=1

extractindexfromfile.pl -in list.notin.YEATS -idx 1 
yeats_MakeNRSet.csh list.notin.YEATS.0 DB/DB.list.notinYEATS  FASTADIR_MAKER notin.YEATS 60 
parseBlastLatestList.pl -out BLASTOUTSUBSET/notin.YEATS.nonredundant -list list.notin.YEATS.0 -blastdir BLASTOUTSUBSET/ -blastcutof $cutoff -forWGS 0 -findcha 0 -isNT=1



##################################################
### Get Annotations for these grouped nonredundant 
##################################################


extractthoseinlist.pl -inf MAPDIRS/maker2PLANTS.60.anno.real -list BLASTOUTSUBSET/notin.YEATS.nonredundant.60.GROUPED -out BLASTOUTSUBSET/notin.YEATS.nonredundant.60.GROUPED.anno
 extractthoseinlist.pl -inf MAPDIRS/finalanno.60.anno.real -list BLASTOUTSUBSET/notin.MAKER-P.nonredundant.60.GROUPED -out BLASTOUTSUBSET/notin.MAKER-P.nonredundant.60.GROUPED.anno


