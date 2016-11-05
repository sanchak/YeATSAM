#!/bin/csh -f

if($#argv != 4  ) then 
  echo "<list> <ksize> <hopping> <tag>"
  exit 
endif 

set list = $1
set ksize = $2
set hopping = $3
set tag = $4


testkmer.pl -genom  DB/DB.$list -out $tag  -lis $list -ksize $ksize -isNT 0 -fasta FASTADIR_ORFBEST -hopping $hopping 
mappingAddCount.pl -in $tag/$tag.0.$ksize.mapall 
groupBasedonCutoff.pl -in $tag/$tag.0.$ksize.mapall.mapped.pw -out $tag/group.$ksize -cutoff 0 
TW $list kmerMergedNonRetroLRR.Self/group.$ksize.allingroup 
cat $tag/group.$ksize.first ofhinAbutnotinB > ! $tag/geneset.$ksize

wc -l $tag/geneset.$ksize
