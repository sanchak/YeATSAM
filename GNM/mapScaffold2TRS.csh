#!/bin/csh -f

if($#argv != 5  ) then 
  echo "Usage : <trs> <scaffold>  <fastadir> <scaffolddir> <buffersize for slicing> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 

set TRS=$1
set SCAFFOLD=$2
set FASTADIR = $3
set SCAFFOLDDIR = $4
set BUFFERSIZE = $5

set BLASTDIR=BLASTOUT_SCAFFOLD

set SCAFFOLDFASTA=$SCAFFOLDDIR/$SCAFFOLD.ALL.1.fasta
set SCAFFOLDFASTACOMP=$SCAFFOLDDIR/$SCAFFOLD.ALL.1.fasta.comp.fasta

set BLASTFILE = $BLASTDIR/$TRS.$SCAFFOLD.blast
set TRSFASTA = $FASTADIR/$TRS.ALL.1.fasta
if(! -e $TRSFASTA || ! -e $SCAFFOLDFASTA) then 
	echo "Either $TRSFASTA  or $SCAFFOLDFASTA does not exist"
	exit 
endif 

mkdir -p $BLASTDIR; mkdir -p INFO; mkdir -p TMP;

myblastcompare2Fastafiles.csh $TRSFASTA $SCAFFOLDFASTA $BLASTFILE N

set REVFILE=TMP/$TRS.$SCAFFOLD.isrev
\rm -rf $REVFILE
checkBLASTtrs2scaffolddirection.pl -inf $BLASTFILE -outf $REVFILE 

if( -e $REVFILE) then 
	echo "======== Reverese ==========="
	if(-e $SCAFFOLDFASTACOMP) then 
	else
	      echo "doing comp"
          writeCompFasta.pl -in $SCAFFOLDFASTA 
	endif

    myblastcompare2Fastafiles.csh $FASTADIR/$TRS.ALL.1.fasta $SCAFFOLDFASTACOMP  $BLASTFILE N
else
	echo "======== SAME ==========="
endif

\rm -rf $REVFILE


## Order is reverse
mapScaffold2TRS.pl -trs $TRS -sname $SCAFFOLD -SCAFFOLDDIR $SCAFFOLDDIR -buffer $BUFFERSIZE

ls SCAFF.INFO/*$TRS*$SCAFFOLD*
