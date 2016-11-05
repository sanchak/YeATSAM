#!/bin/csh -f

if($#argv != 6  ) then 
  echo "ydiffwithfDIR.csh  <fastaALIST> <fastaADIR> <fastaB_DB> <blast> <cutoff> <tag>"
  exit 
endif 


set fastaALIST = $1
set fastaADIR = $2
set fastaDB = $3
set blast = $4
set cutoff = $5
set tag = $6

mkdir -p $tag 

if($blast == blastp) then
    set comm=makeblastdbP.csh
else
	set comm=makeblastdbN.csh
endif

## create DB

egrep -l "\|" $fastaDB
fixFastaNm.pl -in $fastaDB -out $tag/tmp -same 
if( ! -e $fastaDB.uniq.psi) then 
	uniquifyFasta.pl -in $fastaDB
    $comm $fastaDB.uniq
endif 


mkdir -p $tag/BO_TMP
setupCommandsFromListAndSchedule.pl -lis $fastaALIST -out $tag.exec -sleep 10 -inter 1000 -string "$blast -db $fastaDB.uniq -query $fastaADIR/REPLACE.ALL.1.fasta -out $tag/BO_TMP/REPLACE.blast.nt" -force 

echo parseBlastLatestList.pl -out $tag.diff -lis $fastaALIST -blastdir $tag/BO_TMP -blastcutoff $cutoff -forWGS 1 -findcha 0 -isNT=1 > ! $tag/finalcommand.csh 
wc -l $tag/finalcommand.csh


