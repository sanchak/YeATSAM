### this creates the ANNOTATE files - you need to have at least one residue in the file clasp.in
### the rest will be filled in...
createCLASPinput.csh 2MALA clasp.in 4 1
createCLASPinput.csh 1FK0A clasp.in 4 1 ## coincidentally both have Val7 

ann2simpleinput.pl -in ANNOTATE/2MALA.outconf.annotated -out 2MALA.simple.in


\rm ANNOTATE
ln -s ANNOTATE.4 ANNOTATE 

runRefExtractEasilyNamed.csh l.ref l.query 0

## These are the final files
wc -l Extract.l.ref.l.query/SINGLE*/*



printPairwise.pl -out $i.out -c $CONFIGGRP -ra 222 -pr $i -in Extract.l.ref.l.query/1FK0A/2MALA.pdb.out


### NOW to PREMON

## this creates a premon one from the CLASP ...
createPremoninputFromCLASP.csh $i 4 polar
mv 2MALA.4.polar.premon.in ANNOTATE.4/

## while this creates a premon from just residues numbers - creating CLASP in the process
createPremoninput.csh $i clasp.in 4 polar
mv 2MALA.4.polar.premon.in ANNOTATE.4/

## polar is just a tag...
runPremonClasp.csh l.ref l.ref 1 1 4 $SRC/PREMONCONFIGS/aalist $cwd/ANNOTATE.4/ polar
