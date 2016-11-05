#!/bin/csh -f

if($#argv != 4  ) then 
  echo "Usage : <list> <FASTADIR> <out> <<command>"
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "


set list = $1 
set FASTADIR = $2 
set out = $3 
set command = $4 

set exitfromcomm = 0
   foreach i (`cat $list`)
   	    if(! -e $FASTADIR/$i.ALL.1.fasta) then 
			echo "No file $FASTADIR/$i.ALL.1.fasta" 
			set exitfromcomm = 1
		endif
  end

if ($exitfromcomm == 1) then
	unlink $out
	echo "Errors, so exiting"
    #exit ;
endif 
 
if(! -e $out) then 
   newfile.csh $out 
   foreach i (`cat $list`)
   	    if( -e $FASTADIR/$i.ALL.1.fasta) then 
	        cat $FASTADIR/$i.ALL.1.fasta >> $out
	   endif 
   end 
else
	echo "Info: File $out exists"
endif

wc -l $list
nfasta $out

if($command == "") then 
else
    $command $out 
endif 
