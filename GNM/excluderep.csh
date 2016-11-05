#!/bin/csh -f

if($#argv != 1  ) then 
	die 
endif 

set original  = $1
echo $original 

egrep -v -i "methyltransferases|intregase|cdna. isNT|est. isNT|mrna|abinitio|ATP synthase subunit delta|UDP-Glycosyltransferase|Zinc finger|cdna.est.omrna|Warning:Did |Expressed protein|glutathione S-transferase|domain-containing protein|disease resistance|Ring.U-box|Cytochrome|F-box|splicing|transmembrane|Gag protease|chromosome|polyproteint|repeat|ABC transporter|retro|Predicted protein|polymerase|hypothe|LRR|kinase|unknown|unchar|pep:novel|gag-|reverse transcriptase|transpos|Uncharacterized protein" $original > ! $original.clean

sortOnLast.pl -in $original.clean -idx 2 
sortOnLast.pl -in $original.clean -idx 2 -cutoff 70 
sortOnLast.pl -in $original.clean.sort -idx 1 -cutoff 100 -rev

wc -l $original*
