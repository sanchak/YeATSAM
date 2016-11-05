#!/bin/csh -f

if($#argv != 5  ) then 
  echo "Usage :  ydiffbasicFASTA <fastaA>   <fastaB> <comm> <cutoff> <tag>"
  echo "HEre you give two FASTA files, we will split A in single ones, create a DB from B and BLAST"
  echo "Makes sense to have lesser number in A"
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "


## TODO - cehck A has lesser number of fastas
set fastaA = $1
set fastaB = $2
set blast = $3
set cutoff = $4
set tag = $5

mkdir -p $tag 


if(! -e $fastaA.FDIR) then 
    mkdir -p $fastaA.FDIR
    setenv FASTADIR $fastaA.FDIR
    splitFasta.pl -in $fastaA 
endif 
extractFastaNames.pl -fast  $fastaA -out  $fastaA.list


ydiffBasicLIST.csh $fastaA.list $fastaA.FDIR  $fastaB $blast $cutoff $tag
