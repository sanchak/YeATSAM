#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use MyGeom;
use PDB;
use ConfigPDB;
use Math::Geometry ;
use Math::Geometry::Planar;


use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($datadir,$ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 0 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "datadir=s"=>\$datadir ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);
usage( "Need to give a input file name => option -listfile ") if(!defined $listfile);
usage( "Need to give a input file name => option -datadir ") if(!defined $datadir);

my $CNT = 0 ; 
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();

my $PWD = cwd;
my @list= util_read_list_sentences($listfile);
my $list = {};
map { s/\s*//g ; $list->{$_} = 1 ; } @list ;


my $info = {};

foreach my $XX (@list){
my $infile = "$datadir/$XX.txt";
print "Reading $infile\n" if($verbose);
my $ifh = util_read($infile);
my ($snm ,$id, $Division,$Lineage ) ;
my $startedLineage = 0 ;
my $done = {};
while(<$ifh>){
     next if(/^\s*$/);
	 next if(/^\s*#/);
	 last if(/^\s*Genetic Code/);
     chomp ;

	 if(/^\s*Scientific Name/){
		$_ = <$ifh>;
		my @l = split ;
	    $snm = $l[0];
     }
	 if(/^\s*Taxonomy ID:/){
		$_ = <$ifh>;
		my @l = split ;
	    $id = $l[0];
     }
	 if(/^\s*Division/){
		$_ = <$ifh>;
		my @l = split ;
	    $Division = $l[0];
     }
	 if(/^\s*Lineage/){
	    $startedLineage = 1;
        next ;
     }

	if($startedLineage){
		my @l = split ;
	    $Lineage = $Lineage . $l[0];
     }
}
   my $retstr  = "$snm\t$id\t$Division\t$Lineage";
   $info->{$snm} = $retstr ;
   close($ifh);
}

foreach my $k (keys %{$info}){
	my $v = $info->{$k};
	print $ofh "$v \n";
}


chmod 0777, $outfile ;

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
