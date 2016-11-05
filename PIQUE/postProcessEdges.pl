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
use Memory::Usage;
use Algorithm::Combinatorics qw(combinations) ;
my $mu = Memory::Usage->new();
$mu->record('');
use AAConfig;



use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($fastafile,$idx,$infile,$postfix,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "fastafile=s"=>\$fastafile ,
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
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
if(-z $infile){
	print "Info: input file $infile is empty - possibly no chains. Exiting\n";
	exit(1);
}
my $ifh = util_read($infile);

#usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
#usage( "Need to give a protein pdb id -option -protein  ") if(!defined $protein);

my $CNT = 0 ; 
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();
my $aaconfig = new AAConfig("$SRC/aa.config");

# for reading fasta
#my ($fe,$x) =  util_readfasta($infile);
#chomp $x ;
my $PWD = cwd;

system ("mkdir -p WORKDIR/PYMOLCOMMANDS");
system ("mkdir -p WORKDIR/CLASPIN");

my @expected = qw ( HSP90_C1<->RAR1_C1
HSP90_C1<->SGT1A_C2
HSP90_C2<->RAR1_C2
HSP90_C2<->SGT1A_C1
HSP90_C1<->SGT1A_C1
HSP90_C2<->SGT1A_C1
HSP90_C2<->SGT1A_C2
HSP90_C3<->SGT1A_C2
HSP90_C3<->SGT1A_C3
HSP90_C1<->CDC37_C1
);

my $expectedtable = util_make_table(\@expected);

my $ignoretable ;
if(defined $ignorefile){
    ($ignoretable) = util_maketablefromfile($ignorefile);
}


##my @list= util_read_list_sentences($listfile);
#my $list = {};
#map { s/\s*//g ; $list->{$_} = 1 ; } @list ;
#
#foreach my $i (@list){
#}


my $FULLIDS = {};
my $PREFIX_2_FULLIDS = {};
my $ID_2_EDGES = {};
## Read the contact file...
my $PDB ;
while(<$ifh>){
     next if(/^\s*$/);
     next if(/^\s*#/);
     chomp ;
	 if(/EDGE/){
	 	s/=//g;
	    my ($x1,$jjj1,$jjj2,$A,$B,$H1,$H2,$res1,$res2) = util_ParseContactFile($_);
	    my ($x2,$PDB1,$PDB2,$copyA,$copyB) = util_ParseContactFile($_);
		if($copyA =~ /HSP/ && $copyB =~ /HSP/){
			print "$_ \n"; ; 
		}
	 }
}


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
