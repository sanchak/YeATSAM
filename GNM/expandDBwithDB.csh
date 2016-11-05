#!/bin/csh -f

if($#argv != 5  ) then
  echo "Usage : <DB1toexpand> <DB2> <FASTADIR> <kmersize> <tag>"
  exit
endif
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

set DB1 = $1
set DB2 = $2
set FASTADIR_DB2 = $3
set kmersize = $4
set tag = $5

mkdir -p $tag
   
   
if(! -e $tag/$tag.0.$kmersize.mapall) then 
    setenv FASTADIR $FASTADIR_DB2
    uniquifyFasta.pl -in $DB1 -write 1 -out $tag/DB1

    #testkmer.pl -genom $DB2 -out $tag -lis $tag/DB1.list -ksize $kmersize -isNT 0 -fasta $FASTADIR_DB2
    testkmerALLMasks.pl -genom $DB2 -out $tag -lis $tag/DB1.list -ksize $kmersize -isNT 0 -fasta $FASTADIR_DB2
	source $tag/commands.csh
endif 

cat $tag/$tag.*mapall > ! $tag/MAPALL

expandDBwithDBCreateNewDB.pl -in $tag/MAPALL -out $tag/DB.$tag -fast $FASTADIR_DB2/ -orig $tag/DB1.list

nfasta $tag/DB.$tag 

