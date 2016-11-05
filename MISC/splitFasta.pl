#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use MyUtils;


my ($infile,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions);
my $howmany = 100000 ;
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
my $ifh = util_read($infile);
my $info = {};
my $ofh ;
system ("mkdir -p $FASTADIR");
print "Writing to $FASTADIR\n";

my $ofhlen = util_write("$infile.length");

my $STR = "";
my $NAME = "";


my $kmertable = {};
while(<$ifh>){
	 next if(/^\s*$/);
     if(/^\s*>/){

	 	## process one fast
	 	if(defined $ofh){
			my $len = length($STR);
			print $ofhlen "$NAME $len\n";
            $STR = "";
	 	    close($ofh) ;
		}

		my @l = split ;
		$NAME =  $l[0];
		$NAME =~ s/>//;
        $ofh = util_write("$FASTADIR/$NAME.ALL.1.fasta");
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




system("numfiles $FASTADIR");

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
