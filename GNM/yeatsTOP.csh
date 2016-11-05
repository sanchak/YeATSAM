#!/bin/csh 

if($#argv != 3  ) then 
  echo "Usage :yeatsTOP.csh <fastafile> "
  exit 
endif 

set initialfile = $1
set genome=$2
set TRS=$3


#dos2unix $1 
mkdir -p FASTADIR_NT/
mkdir -p ORF/
mkdir -p ORFTMP/
mkdir -p FASTADIR_ORF/

#### DB/all.fasta is the initial files - make sure you have the name correct for the TRS's, 
#### First get a list of trs - remove exact ones.
#### for example, no semicolons, etc
set initialfileuniq = $initialfile.uniq 
set initialfileuniqsam = $initialfileuniq.sam
\cp -f $1 $initialfile





setenv FASTADIR $cwd/FASTADIR_NT/

if(! -e $initialfileuniq) then 
   uniquifyFasta.pl -in $initialfile -out $initialfileuniq -write 1 
endif 


orfprocessList.csh $initialfileuniq.list



if(! -e $genome.bwt) then 
   bwa index $genome 
endif 
if(! -e $initialfileuniqsam.unmapped.genome.sam) then 
    bwa mem $genome $initialfileuniq > ! $initialfileuniqsam
    samtools view -f4 -S $initialfileuniqsam  > ! $initialfileuniqsam.unmapped.genome.sam
    samtools view -F4 -S $initialfileuniqsam  > ! $initialfileuniqsam.mapped.genome.sam
else
	echo " $initialfileuniqsam.unmapped.genome.sam exists"
endif


if(! -e $TRS.bwt) then 
   bwa index $TRS 
endif 
if(! -e $initialfileuniqsam.unmapped.trs.sam) then 
    bwa mem $TRS $initialfileuniq > ! $initialfileuniqsam
    samtools view -f4 -S $initialfileuniqsam  > ! $initialfileuniqsam.unmapped.trs.sam
    samtools view -F4 -S $initialfileuniqsam  > ! $initialfileuniqsam.mapped.trs.sam
else
	echo " $initialfileuniqsam.unmapped.trs.sam exists"
endif 





