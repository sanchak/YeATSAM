#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use MyGeom;
use MyGNM;
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
my ($fastafile,$idx,$annotatefile,$postfix,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "fastafile=s"=>\$fastafile ,
            "postfix=s"=>\$postfix ,
            "annotatefile=s"=>\$annotatefile ,
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
usage( "Need to give a input file name => option -annotatefile ") if(!defined $annotatefile);
my $ifh = util_read($annotatefile);
usage( "Need to give a input file name => option -fastafile ") if(!defined $fastafile);
my $ifhfasta = util_read($fastafile);

$outfile = "$fastafile.anno" if(!defined $outfile);
my $ofh = util_write($outfile);

my $info = {};
while(<$ifh>){
     next if(/^\s*$/);
	 next if(/^\s*#/);
	 my ($nm) = split ; 
	 $info->{$nm} = $_ ;
}
close($ifh);

while(<$ifhfasta>){
     next if(/^\s*$/);
	 next if(/^\s*#/);
	 if(/^>/){
	 	s/>//;
		my ($nm) = split ;
		if(exists $info->{$nm}){
			my $ann = $info->{$nm} ;
			$ann =~ s/#none none//;
			print $ofh ">$ann";
		}
		else{
			print $ofh ">$nm\n";
		}
	 }
	 else{
	 	print $ofh $_ ;
	 }
}


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
