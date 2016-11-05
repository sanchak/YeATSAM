#!/bin/csh -f

if($#argv != 3  ) then 
  echo "Usage : pdb hetatm ispolar"
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 

set pdb=$1
set hetatm=$2
set ispolar=$3


hetatmProc.pl -outf ooo -l l -het $hetatm -ispolar $ispolar -p $pdb

