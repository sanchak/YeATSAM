#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use MyUtils;


my ($infile,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions);
my $howmany ;
my $verbose = 1 ;
GetOptions(
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "which_tech=s"=>\$which_tech ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);

my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a input file name => option -howmany ") if(!defined $howmany);
my $ifh = util_read($infile);
my $info = {};
my $ofh ;


my ($llll,$nfasta) = `nfasta $infile ` ;
my $nperfile = int($nfasta / $howmany) ;


print "Divide $infile with $nfasta into $howmany different files with $nperfile \n";
my $STR = "";
my $NAME = "";
foreach my $i (0..$howmany){
         util_write("$infile.$i");
}


my $kmertable = {};
my $CNT = 0 ;
my $counter = 0 ;
while(<$ifh>){
	 next if(/^\s*$/);
     if(/^\s*>/){

	 	## process one fast
		$CNT++ ;
		if($CNT eq $nperfile){
			$CNT = 0 ;
			$counter++;
		}
		my @l = split ;
		$NAME =  $l[0];
		$NAME =~ s/>//;
        $ofh = util_open_or_append("$infile.$counter");
		print $ofh $_ ;
	 }
	 else{
	    die if(!defined $ofh);
	    die if(/>/);
		## remove spaces 
		s/ //g;
	 	print $ofh $_ ;
		chomp;
		$STR = $STR . $_ ;
	 }
}
close($ifh);



sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
