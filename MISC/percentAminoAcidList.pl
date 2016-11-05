#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use MyGeom;
use PDB;
use ConfigPDB;

use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($commands,$percentfile,$infile,$outfile,$template,$listfile,$protein);
my (@expressions,$config);
my $aminoacid = 1;
my $force = 0;
my $verbose = 1 ;
my $aadiff = "aadiff.rep";
GetOptions(
            "template=s"=>\$template ,
            "commands=s"=>\$commands ,
            "percentfile=s"=>\$percentfile ,
            "aadiff=s"=>\$aadiff ,
            "protein=s"=>\$protein ,
            "config=s"=>\$config ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "aminoacid=i"=>\$aminoacid ,
            "force=i"=>\$force ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_open_or_append($outfile);
usage( "Need to give a listfile pdb id -option -listfile  ") if(!defined $listfile);
usage( "Need to give a config file name => option -config ") if(!defined $config);
usage( "Need to give a commands file name => option -commands ") if(!defined $commands);

my $CNT = 0 ; 
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my $PWD = cwd;


my ($donetable) = util_maketablefromfile($outfile);
my @list= util_read_list_sentences($listfile);
my $ofhcommands = util_write($commands);
foreach my $protein (@list){
    if(! exists $donetable->{$protein} && !$force){
	    print $ofhcommands "percentAminoAcid.pl -pro $protein -out $outfile -inf $FASTADIR/$protein.ALL.1.fasta -con $config\n";
    }
}



sub usage{
    my ($msg) = @_ ;
	    print $msg , "\n" ;
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
		    die ;
			}
