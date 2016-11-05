#!/bin/csh -f

if($#argv != 1  ) then 
	echo "hossein.csh <genome>"
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "


set listI = INPUTFILES/list.I
set listN = INPUTFILES/list.N

newfile.csh BOTHlists
cat $listI $listN >> BOTHlists

set genome = $1

## INDEX the genome
if(! -e $genome.amb) then 
    echo bwa index $genome
endif 

mkdir -p SAMFILES
mkdir -p SAMFILESPW

### EACH BLOCK CAN BE RUN in parallel

## MAP to genome
foreach i (`cat BOTHlists`)
  if(! -e SAMFILES/$i.sam) then
     echo "bwa mem $genome  INPUTFILES/$i.R1.fa INPUTFILES/$i.R1.fa   > ! SAMFILES/$i.sam "
  endif 
end 

## Diff out non-genomic fragements
foreach i (`cat BOTHlists`)
	set what = unmapped 
    if(! -e SAMFILES/$i.$what.sam) then 
       hosseinSingleSam.csh $i SAMFILES $what > & ! /dev/null & 
	endif 
end 



newfile.csh PW.csh
foreach i (`cat $listI`)
     foreach j (`cat $listN`)
	 	if(! -e SAMFILESPW/$i.$j.sam ) then 
          echo "bwa mem SAMFILES/$i.unmapped.fa  SAMFILES/$j.unmapped.fa  > ! SAMFILESPW/$i.$j.sam " >> PW.csh 
	    endif 
	 end

	 # only one infected
	 #break
end 

foreach i (`cat $listI`)
     foreach j (`cat $listN`)
	 	if(! -e SAMFILESPW/$i.$j.unmapped.sam ) then 
          hosseinSingleSam.csh $i.$j SAMFILESPW  unmapped & 
          hosseinSingleSam.csh $i.$j SAMFILESPW  mapped  & 
	    endif 
	 end

	 # only one infected
	 break ;
end 

