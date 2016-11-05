#!/bin/csh -f
if($#argv != 1  ) then 
  echo "Usage <TAG>"
  exit 
endif 

set TAG = $1

## PLEASE CHECK FOR LOWERCASE, symbols, etc...
foreach i ( BLASTOUT_VIRAL BLASTOUT_AI OUTDIR  BLASTOUT_2WGS  FASTADIR_NT  BLASTOUT_AA  BLASTOUT_2WGS  LISTS BLASTOUT_RRNA FASTADIR_ORFBEST.STARTCODON )
	mkdir -p $i
end 


## PLEASE CHECK FOR LOWERCASE, symbols, etc in the TRSFILE

#set TRSFILE = DB/$TAG.trs.fasta
#set TRSFILE = FASTADIR_NT/$TAG.ALL.1.fasta

set LIST =  $TAG.list
set LISTCLEAN =  $LIST.clean
set ORFS3 = $LIST.orfall3
set LISTCLEANBLAST =  $LIST.clean.BLAST
set ANNFILE =  $LISTCLEANBLAST.anno #TODO - extract
set LISTCLEANBLASTANNO =  $ANNFILE.0 #TODO - extract
set LISTCLEANBLASTANNOMERGED  =  $LISTCLEANBLASTANNO.merged
set LISTCLEANBLASTANNOMERGEDUNIQ =  $LISTCLEANBLASTANNOMERGED.uniq

set LISTCLEANBLASTNOANNO =  $LISTCLEANBLAST.noanno

#set DB_WGS = DB/wgs.5d.scafSeq200+.trimmed 
set DB_WGS = DB/chickpea.genome.fasta 
set DB_TAIR = DB/TAIR10_pep_20101214

set DB_ORFS3 = DB/DB.$ORFS3
set DB_ANNO = DB/DB.$LISTCLEANBLASTANNO
set DB_ANNOuniqMERGED = DB/DB.$LISTCLEANBLASTANNOMERGEDUNIQ
set DB_ALL_BEST = DB/DB.$LISTCLEANBLAST.list
set DB_NOANNO = DB/DB.$LISTCLEANBLASTNOANNO

set DB_RRNADB = DB/rel_con_pln_r124.rRNA.fasta 
set DB_MAKERDB   = DB/walnut.wgs.5d.all.maker.proteins.fasta
set DB_REFPROTEOME = DB/plantpep.fasta 

set BO_AA = BLASTOUT_AA





#set OUTPERCENTGOOD=$INITLIST.signalp.$SIZE.annovalues.sort.PERCENTGOOD

set FA_AI = FASTADIR_ORFBEST.STARTCODON
set CUTOFF_AI_SIZE = 60
set SIZESIGNNALP=20
set BO_AI = BLASTOUT_AI
set LISTCLEANBLASTNOANNO_AI = $LISTCLEANBLASTNOANNO.info.sort.0
set DB_AI = DB/DB.$LISTCLEANBLASTNOANNO_AI

