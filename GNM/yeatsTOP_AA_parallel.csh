#!/bin/csh 


#### DB/all.fasta is the initial files - make sure you have the name correct for the TRS's, 
#### First get a list of trs - remove exact ones.
#### for example, no semicolons, etc

set initialfile = DB/all.fasta

mkdir -p FASTADIR_ORFLONGEST
mkdir -p FASTADIR
mkdir -p ORF/
mkdir -p ORFTMP/
mkdir -p FASTADIR_ORF/


# 1) Get list of TRSs, removing exact ones
setenv FASTADIR $cwd/FASTADIR
newfile.csh mapname2length
pruneSameSequenceFromMadeFasta.pl -outf list.trs -inf DB/all.fasta -writedata 1 -length mapname2length
extractindexfromfile.pl -in list.trs -out ttt -idx 0 -sep ""
\mv -f ttt list.trs

## list.trs.mapping.sort has the exact same ones
sort.pl -in list.trs.mapping -idx 1 > & ! /dev/null 



newfile.csh RRR
foreach i (`cat list.unique`)
      echo "getorf FASTADIR/$i.1.ALL.1.fasta ORFTMP/$i.1.orf  > & \! /dev/null & " >> RRR 
end
scheduleprocessInsertingsleep.pl -inf RRR -sleep 1 -inter 500

newfile.csh SSS
foreach i (`cat list.unique`)
   echo "replacestring.pl -wit ".ORF_" -whic "_" -inf  ORFTMP/$i.orf -outf ORF/$i.orf   > & \! /dev/null & " >>SSS 
end
scheduleprocessInsertingsleep.pl -inf SSS -sleep 1 -inter 100


newfile.csh TTT
newfile.csh list.orflongest
foreach i (`cat list.unique`)
   echo "findrepeatedorfs.pl -trs $i -orfdir ORF/ -write 4  > & \! /dev/null & " >>TTT 
end
scheduleprocessInsertingsleep.pl -inf TTT -sleep 1 -inter 100









