
my $SELF = "trinity.list.clean.BLAST.anno.0.merged.uniq.mfix.substring.nonretroLRR";
foreach my $i (3..40){
      my $c = $i * 2 ;
      print "squishOneKmer.csh $SELF $c 0 kmerMergedNonRetroLRR.Self\n";
}
