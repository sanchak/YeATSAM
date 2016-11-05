

source $SRC/GNM/step.setup.csh 



TW $LISTCLEANBLAST.list $ANNFILE
\mv -f ofhinAbutnotinB $LISTCLEANBLAST.noanno


extractindexfromfile.pl -in $ANNFILE
\rm -f $DB_ANNO
createDBFasta.csh $LISTCLEANBLASTANNO  FASTADIR_ORFBEST  $DB_ANNO  makeblastdbP.csh 
uniquifyFasta.pl -in $DB_ANNO

