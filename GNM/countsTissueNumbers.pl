#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use MyGeom;
use MyConfigs;
use PDB;
use ConfigPDB;
use Math::Geometry ;
use Math::Geometry::Planar;


use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($idx,$infile,$postfix,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
my $start = 5 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "postfix=s"=>\$postfix ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "idx=i"=>\$idx ,
            "verbose=i"=>\$verbose ,
            "start=i"=>\$start ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a input file name => option -start ") if(!defined $start);
my $ifh = util_read($infile);

#usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
#usage( "Need to give a protein pdb id -option -protein  ") if(!defined $protein);
my $CNT = 0 ; 
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();

# for reading fasta
#my ($fe,$x) =  util_readfasta($infile);
#chomp $x ;
my $PWD = cwd;


print "start = $start, this defines whether its anno\n";
my $ignoretable ;
if(defined $ignorefile){
    ($ignoretable) = util_maketablefromfile($ignorefile);
}

my $sortNumberofTissues = {};
my $sortMax = {};

my $finalN = 0 ;
while(<$ifh>){
	 my ($nm,@l) = split ; 
	 $finalN++;
	 foreach my $idx (1..$start){
	 	shift @l ;
	 }
	 my $N = @l ;
	 my $max = util_get_max(\@l);
	 my $tenpercent = ($max / 10) - 50 ; ## the 5 is to get rid of smalleer numbers

	 my $cnt = 0 ;
	 foreach my $i (@l){
	 	if($i && $i > 50 && $i >= $tenpercent){
			$cnt++;
		}
	 }
	 #print "$nm $max $cnt\n";
	 $sortNumberofTissues->{$_} = $cnt ;
	 $sortMax->{$_} = $max ;
}

my @s = sort { $sortNumberofTissues->{$a} <=> $sortNumberofTissues->{$b}} ( sort keys %{$sortNumberofTissues});

print $ofh " TRS   CE   CI   CK   EM   FL   HC   HL   HP   HU   IF   LE   LM   LY   PK   PL   PT   RT   SE   TZ   VB\n";
if(defined $cutoff){
    my @next ;
    foreach my $i (1..$cutoff){
	    my $x = shift @s ;
	    push @next, $x ;
    }
    
    my @sortednext = sort { $sortMax->{$b} <=> $sortMax->{$a}} (@next);
    
    print $ofh @sortednext;
}
else{
	print $ofh "@s";
}


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
