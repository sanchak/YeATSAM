#!/bin/csh -f

if($#argv != 1  ) then 
  #echo "<list> <outfile> <writeind>"
  echo "<list> "
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

set list = $1 
#set outfile = $2 
#set writeind = $3 
#touch $outfile

echo DSSP=$DSSP and HELIXDIR=$HELIXDIR
touch HTH
foreach i (`cat $list`)
	if(! -e "$DSSP/$i.dssp" || -z "$DSSP/$i.dssp") then
	    echo mkdssp -i $PDBDIR/$i.pdb -o $DSSP/$i.dssp
	    mkdssp -i $PDBDIR/$i.pdb -o $DSSP/$i.dssp
	else
		echo Info: DSSP $DSSP/$i.dssp file exists
	endif


	## Do not process individual files here ...

	#grep -w $i $outfile > & ! /dev/null 
	#if  ($status != 0 ) then 
	   #helixparseDSSPoutput.pl -outf $outfile -p $i -dssp $DSSP/$i.dssp -what BETA -writeind $writeind
	   #helixparseDSSPoutput.pl -outf $outfile -p $i -dssp $DSSP/$i.dssp -what HELIX -writeind $writeind
	#else
		#echo $i already processed in output file $outfile
	#endif 

end 

 

