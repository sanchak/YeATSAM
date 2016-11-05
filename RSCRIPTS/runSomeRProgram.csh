#!/bin/csh  -f

if($#argv != 5  ) then 
    echo "Wrong args, required 4 " ; exit 
endif 

set infile=$1
set idx1=$2
set idx2=$3
set printcom=$4
set meantotest=$5
set outfile=$infile.commands


newfile.csh $outfile
runSomeRProgram.pl -inf $infile -var A -outf $outfile -idx $idx1 -scal 1 -mean $meantotest
runSomeRProgram.pl -inf $infile -var B -outf $outfile -idx $idx2 -scal 1 -printcom  $printcom -mean $meantotest

R --no-save < $outfile  > & ! $outfile.log
