#!/usr/bin/perl -w 
my $postfix = "pdb";
my @l = <*.$postfix>;
my $N = @l ;
print "Found $N files with postfix $postfix\n";

foreach my $i (@l){
  my $orig = $i ;
  $i =~ s/\s*//g;
  $i = uc($i);
  $i =~ s/\.PDB/\.pdb/;
  if($orig ne $i){
      system ("mv \"$orig\" $i");
  }
  else{
  	print "no spaces in $i\n";
  }
}
