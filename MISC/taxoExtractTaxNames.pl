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
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();

my $info = {};
my $done = {};

my $ofhignored = util_write("ignored");
my $ofhlist = util_write("list.meta");
while(<$ifh>){
	my $orig = $_ ;
     next if(/^\s*$/);
	 next if(/^\s*#/);
	 s/PREDICTED://;
	 s/TPA_inf://;
	 if(/Uncultured/i || /E.coli/ || /Escherichia/ || /E. coli/ || /Bacterium/ || /Single read/ || /Synthetic construct/ || /Uncultured eukaryote/i || /Cloning vector/i){
	 	print $ofhignored "$_";
		next ;
	 }
	 my ($x,$y,$a,$b) = split ;

	 if($a =~ /\./){
	 	 my $lll = $a;
	 	 $lll =~ s/\./  /;
	     ($a,$b) = split " ",$lll ;
	     my $str = "$b";
	     print $ofhlist "$orig";
		 next if(exists $done->{$str});
		 $done->{$a} = 1 ;
		 $done->{$b} = 1 ;
	     print $ofh "$str\n";
	 }
	 else{
	 	 $b = "" if($b eq "sp.");
	     my $str = "$a $b";
	     print $ofhlist "$orig";
		 next if(exists $done->{$a});
		 $done->{$a} = 1 ;
		 $done->{$b} = 1 ;
	     print $ofh "$a\n";
	     print $ofhlist "$orig";
	 }

}

system("wc -l ignored list.meta");

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
