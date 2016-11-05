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
my $ofh = util_write("$infile.stats") ;
system ("mkdir -p $FASTADIR");


my $STR = "";
my $NAME = "";

my $nA = 0 ; 
my $nT = 0 ; 
my $nG = 0 ; 
my $nC = 0 ; 

my $CNTfasta = 0 ;
while(<$ifh>){
	 next if(/^\s*$/);
     if(/^\s*>/){
	 	 $CNTfasta++;
	 }
	 else{
		s/ //g;
		chomp;

	    my @A = (/(A)/g);
	    my @T = (/(T)/g);
	    my @G = (/(G)/g);
	    my @C = (/(C)/g);
		$nA = @A + $nA ;
		$nT = @T + $nT ;
		$nG = @G + $nG  ;
		$nC = @C + $nC ;
	 }
}

my $sum = $nA + $nT + $nG + $nC ;
my $sumAT = $nA + $nT ;
my $sumGC = $nG + $nC ;
my $nATPercent = int (($sumAT/$sum) * 100) ;
my $nGCPercent = 100 - $nATPercent ;
print $ofh "$infile CNTfasta = $CNTfasta. sum = $sum, $nA + $nT + $nG + $nC \n" ;
print $ofh "$infile AT= $nATPercent%, GC = $nGCPercent \n";

close($ifh);





sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
