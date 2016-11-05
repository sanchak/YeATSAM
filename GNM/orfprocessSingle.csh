#!/bin/csh -f
if($#argv != 2  ) then 
	echo "usage : <i> <FASTADIR>"
	exit 
endif 
mkdir -p ORF
mkdir -p ORFTMP
set i=$1 
set FDIR=$2

if(! -e ORF/$i.orf) then 
       getorf $FDIR/$i.ALL.1.fasta ORFTMP/$i.orf 
       fixORFnames.pl -inf ORFTMP/$i.orf -out ORF/$i.orf
endif 
findrepeatedorfs.pl -trs $i -orfdir ORF/ -write 4 

