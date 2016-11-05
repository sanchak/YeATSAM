#!/bin/csh -f

if($#argv != 2  ) then 
  echo "Usage: <trs> <scaffold>"
  exit 
endif 

set TRS=$1
set SCAFFOLD=$2

mkdir -p BLASTOUT; mkdir -p LOGS; mkdir -p NEWFASTA; mkdir -p SCAFFOLDDIR; mkdir -p INTRONS

touch infointrons.csv 

setenv FASTADIR_TRS /home/sandeepc/DATA/rafael/walnut/2/FASTADIR_ORFNEW/
setenv FASTADIR /home/sandeepc/DATA/rafael/walnut/2/NEW/FASTADIR_NT
setenv BLASTDB /home/sandeepc/DATA/rafael/walnut/DB

set    GENOMEDB=wgs.5d.scafSeq200+.trimmed
set  SCAFFOLDDIR=SCAFFOLDDIR



set SCAFFOLDFASTA=$SCAFFOLDDIR/$SCAFFOLD.ALL.1.fasta
setenv FORCERUN 1
#unsetenv FORCERUN 
if( $?FORCERUN  ||  ! -e $SCAFFOLDDIR/$SCAFFOLD.ALL.1.fasta) then 
    extractsinglefastafromfile.pl -inf $BLASTDB/$GENOMEDB -wh $SCAFFOLD
else
	echo "$TRS and $SCAFFOLD seems to be processed. Remove $SCAFFOLDDIR/$SCAFFOLD.ALL.1.fasta  to force"
	#exit 
endif



\rm -rf $TRS.isrev
myblastcompare2Fastafiles.csh $FASTADIR/$TRS.ALL.1.fasta $SCAFFOLD.ALL.1.fasta BLASTOUT/$TRS.blast N


parseBlastLatestpairwise.pl -outf lll -inf BLASTOUT/$TRS.blast -trs $TRS -ver 1 -justchecking -scaff $SCAFFOLDFASTA

if( -e $TRS.isrev) then 
	echo "Doing the reverse"
    writeCompFasta.pl -in $SCAFFOLD.ALL.1.fasta
	\mv -f $SCAFFOLD.ALL.1.fasta.comp.fasta $SCAFFOLDDIR/$SCAFFOLD.ALL.1.fasta 
	unlink $SCAFFOLD.ALL.1.fasta
else
	\mv -f $SCAFFOLD.ALL.1.fasta $SCAFFOLDDIR/$SCAFFOLD.ALL.1.fasta 
endif


### This is the real thing - no justchecking....
myblastcompare2Fastafiles.csh $FASTADIR/$TRS.ALL.1.fasta $SCAFFOLDDIR/$SCAFFOLD.ALL.1.fasta BLASTOUT/$TRS.blast N
parseBlastLatestpairwise.pl -outf $TRS.$SCAFFOLD.info -inf BLASTOUT/$TRS.blast -trs $TRS -ver 1 -scaff $SCAFFOLDFASTA


\rm -rf $TRS.isrev
