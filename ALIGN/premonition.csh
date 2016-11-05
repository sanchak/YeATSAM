#!/bin/csh -f

if($#argv != 3  ) then 
  echo "Usage : <list> <mlen> <config>"
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 

set PWD = ` pwd`
set WORKDIR = $PREMONITION
set list = $PWD/$1
set mlen = $2
set config = $PWD/$3

echo WORKDIR is $WORKDIR
mkdir -p $WORKDIR/$mlen
cd $WORKDIR/$mlen
foreach i (`cat $list`)
	if(! -e "$i.premon.out.zip") then 
	   # do not parallelize this the per residue - rather parallize the per protein...
       $SRC/ALIGN/premonition.pl -p1 $i  -ml $mlen -list $config 
	   $SRC/ALIGN/premonitionMerge.pl -out $i.premon.out -p1 $i
	   if(! -z $i.premon.out) then 
	      zip -r $i.premon.out.zip $i.premon.out
	   endif 
	   \rm $i*singleout
	   \rm $i.premon.out
   endif 
end 


