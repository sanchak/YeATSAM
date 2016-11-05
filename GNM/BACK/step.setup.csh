#!/bin/csh -f

## PLEASE CHECK FOR LOWERCASE, symbols, etc...
foreach i (  BLASTOUT_2WGS  FASTADIR_NT  BLASTOUT_AA  BLASTOUT_2WGS  LISTS BLASTOUT_RRNA )
	mkdir -p $i
end 

## PLEASE CHECK FOR LOWERCASE, symbols, etc in the TRSFILE
set TRSFILE = DB/Combined_TrinityFull.fasta
set WGS = DB/wgs.5d.scafSeq200+.trimmed 
set RRNADB = DB/rel_con_pln_r124.rRNA.fasta 

set TAG = trinity
set LIST =  $TAG.list






