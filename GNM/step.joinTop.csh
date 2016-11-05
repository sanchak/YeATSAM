s $SRC/GNM/step.join.csh MERGEDIR trinity.list.clean.BLAST.anno.uniq DB/DB.trinity.list.clean.BLAST.anno.0.uniq
createDBFasta.csh MERGEDIR/newlist.1 FASTADIR_ORFBEST/ MERGEDIR/newDB junkcomm

s $SRC/GNM/step.join.csh MERGEDIR.1 MERGEDIR/newannfile.1 MERGEDIR/newDB
