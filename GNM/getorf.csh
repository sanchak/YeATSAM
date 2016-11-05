#!/bin/csh -f
set i=$1 

if(! -e ORF/$i.orf) then 
       getorf $FASTADIR/$i.ALL.1.fasta ORFTMP/$i.orf 
       fixORFnames.pl -inf ORFTMP/$i.orf -out ORF/$i.orf
endif 
