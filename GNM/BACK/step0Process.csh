#!/bin/csh -f

source $SRC/GNM/step.setup.csh 

set tag=trs2scaffold
set blastdir=BLASTOUT_2WGS
set blastcutoff=75
set list = trinity.list
set expectedout = "LISTS/$tag/$tag.$blastcutoff.list"
if( ! -e $expectedout) then
   parseBlastLatestList.pl -out $tag -lis $list -blastdir $blastdir -blastcutoff $blastcutoff -forWGS 1 -findcha 0 -isNT=1  
endif 
## after running we might not get the output, so it should be a diff if
if( ! -e $expectedout) then
    echo "File expectedout: $expectedout not made, quitting"
	exit 
endif 
\cp -f $expectedout trinity.list.clean 

set tag=rrna
set blastdir=BLASTOUT_RRNA/
set blastcutoff=1000
set list = $expectedout
set expectedout = "LISTS/$tag/$tag.$blastcutoff.list"
## after running we might not get the output, so it should be a diff if
if( ! -e $expectedout) then
   parseBlastLatestList.pl -out $tag -lis $list -blastdir $blastdir -blastcutoff $blastcutoff -forWGS 1 -findcha 0 -isNT=1  
endif 
if( ! -e $expectedout) then
    echo "File expectedout: $expectedout not made, quitting"
	exit 
endif 


twolists.pl trinity.list.clean $expectedout
\mv -f ofhinAbutnotinB trinity.list.clean 




