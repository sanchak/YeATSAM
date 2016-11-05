## create the list of 5 lettered PDB - includes the chain

## then make the DSSP files
$SRC/HELIX/dssplist.csh list5

## you can write the individual files thus:
set i=3BIPA
helixparseDSSPoutput.pl -outf fdffsd -p $i -dssp $DSSP/$i.dssp -what BETA -writeind 1
helixparseDSSPoutput.pl -outf fdffsd -p $i -dssp $DSSP/$i.dssp -what HELIX -writeind 1

# In order to process one helix, you can do this...
set i=3CPWW.HELIX2
setenv PDBDIR $HELIXDIR/
helixparseDSSPoutput.pl -outf fdffsd -p $i -dssp $DSSP/$i.dssp -what HELIX -writeind 1
## this would have created a "in" file in tex...
cd FigHelix/

##%%%%@@@@@@@@@@@@##
#There is a problem of handling single helices out of the protein...sometimes DSSP gives different helices...
#So, just use the sequence

set i=3CPWW.HELIX2
setenv PDBDIR $HELIXDIR/
setenv FASTADIR $cwd
writeFasta.pl -p $i 
helixwheel.pl -pro $i -justs $FASTADIR/$i.ALL.1.fasta


## This writes all the info for the helices and beta sheets - ignoring already processed files
processOnePDBforallPropertieslist.csh list5 finaloutfile


## this parallizes in each file
setupCommandsFromListAndSchedule.pl -lis tmp -out fff -sleep 1 -inter 100 -string '$SRC/DB/processOnePDBforallProperties.pl -clean 1 -pro $i -out $HELIXDIR/$i.helix'
