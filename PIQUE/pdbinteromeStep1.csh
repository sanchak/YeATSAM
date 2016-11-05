#!/bin/csh -f

if($#argv != 1  ) then 
  echo "Usage : ./generaterep.csh  <impl dir> <file_having_list_of_designs> <tech - eg altera> <mode> <dirfortech - eg stratixii> "
  echo "You need to set ,  BENCH_HOME , BIN_HOME &  MGC_HOME "
  exit 
endif 
#echo "PERLLIB = $PERLLIB ,  $BENCH_HOME = BENCH_HOME , BIN_HOME = $BIN_HOME , MGC_HOME = $MGC_HOME "

setenv BLASTDB /media/sandeepc/USB30FD/PDBAA/

mkdir -p WORKDIR

set listquery=$1

newfile.csh listofPDB
newfile.csh tmplist
foreach i  (` cat $listquery `)
   #unlink WORKDIR/PDBRESULTS.$i
   if(! -e WORKDIR/PDBRESULTS.$i) then 
       BLASTP $BLASTDB/pdbaa $FASTADIR/$i.ALL.1.fasta WORKDIR/PDBRESULTS.$i
   endif 
   extractHomologousPDBS.pl -outf WORKDIR/listofPDB.$i -in WORKDIR/PDBRESULTS.$i -cutoff 100
   cat WORKDIR/listofPDB.$i >> listofPDB 
   echo "WORKDIR/listofPDB.$i.withchain" >> tmplist 
end


getCommonBetweenFiles.pl -lis tmplist -out ignore -hack 2LY8


newfile.csh listofPDB5
UNIQ listofPDB 
foreach i  (` cat listofPDB.uniq `)
    pdbgetlistAndSplit.csh $i
	cat $PDBDIR/$i.list >> listofPDB5
end

