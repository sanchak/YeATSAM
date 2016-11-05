#!/bin/csh 

if($#argv != 2  ) then 
  echo "Usage: DB mergedir"
  exit 
endif 

## all files are created within the MERGED dir
set db=$1 
set mergedir=$2 

mkdir -p $mergedir

newfile.csh alreadydone
newfile.csh $mergedir/joinscript.1
newfile.csh $mergedir/joinscript.2
foreach i ( 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19  )
     echo  "joinTRS.pl -inf $db -outf data.MERGED.OVERLAP -size $i -mergedir $mergedir & " >> $mergedir/joinscript.1
end


echo "\\rm merge.list alreadydone mapping.merge merge.include merge.exclude " >> $mergedir/joinscript.2
foreach i ( 5 6 7 8 9 10 11 12 13 14 15  16 17 18 19 )
      #joinTRSPostProcess.pl -inf data.MERGED.OVERLAP -size $i -anno annotations.txt
      echo joinTRSPostProcess.pl -fastadir FASTADIR_ORFBEST -inf data.MERGED.OVERLAP -size $i  -mergedir $mergedir -anno map.TRS2scaffold >> $mergedir/joinscript.2
end


