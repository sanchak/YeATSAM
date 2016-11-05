#!/bin/csh -f

if($#argv != 3  ) then 
	echo "hosseinSingleSam.csh <nm>"
  exit 
endif 

set i=$1 
set dir=$2 
set what=$3 

if($what == "unmapped") then 
	set option = f4
else
	set option = F4
endif 

if(! -e $dir/$i.$what.sam) then
     samtools view -$option -S $dir/$i.sam > ! $dir/$i.$what.sam
     sam2fasta.pl -in $dir/$i.$what.sam -out $dir/$i.$what.fa
     bwa index $dir/$i.$what.fa 

else
	echo "Already done $dir/$i.$what.sam ... so doing nothing"
endif 





