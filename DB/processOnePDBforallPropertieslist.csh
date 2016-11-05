#!/bin/csh -f

if($#argv != 2  ) then 
  echo "<list> <outfile>"
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

set list = $1 
set outfile = $2 
touch $outfile

twolists.pl $list $outfile


foreach i (`cat ofhinAbutnotinB`)
	if( -e "$DSSP/$i.dssp") then
       processOnePDBforallProperties.pl -pr $i -out $outfile
	endif 
end 

 

