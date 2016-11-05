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
my ($fastafile,$size,$infile,$postfix,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
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
            "size=i"=>\$size ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a input file name => option -cutoff ") if(!defined $cutoff);
#my ($tablefasta,$infoSeq2PDB,$mapChainedName2Name) = util_parsePDBSEQRESNEW($fastafile,0,0);

#usage( "Need to give a size pdb id -option -size  ") if(!defined $size);

my $PWD = cwd;

my $ifh = util_read($infile);

my $info = {};
my $allinfo = {};
my $fullline = {};
while(<$ifh>){
     next if(/^\s*$/);
	 next if(/^\s*#/);

	 my ($nm,$len,@l) = split ; 
	 $allinfo->{$nm} = 1 ;
	 $fullline->{$nm} = $_ ;
	 my $X = HOWMANYGREATER($nm,$len,15,@l);
	 $info->{$nm} = $_  if($X);

	 $X = HOWMANYGREATER($nm,$len,14,@l);
	 $info->{$nm} = $_  if($X >1);


	 $X = HOWMANYGREATER($nm,$len,12,@l);
	 $info->{$nm} = $_  if($X >2);

	 $X = HOWMANYGREATER($nm,$len,10,@l);
	 $info->{$nm} = $_  if($X >3);

	 $X = HOWMANYGREATER($nm,$len,8,@l);
	 $info->{$nm} = $_  if($X >6);
}

my $NN = (keys %{$info});

my $GOOD = 0 ;
foreach my $k (keys %{$allinfo}){
	if(! exists $info->{$k}){
		my $v = $fullline->{$k};
		print $ofh "$v";
		$GOOD++;
	}
}
print "There were $NN which has percenttage issues, while $GOOD were good\n";


sub HOWMANYGREATER{
	 my ($nm,$len,$CUTOFFPERCENT,@ll) =  @_ ;

	 my $X = 0 ;
	 while(@ll){
	 	my $aa = shift @ll ;
	 	my $percent = shift @ll ;
		next if($aa =~ /L/);
		if($percent > $CUTOFFPERCENT){
			 $X++;
		}
	 }
	 return $X ;
}

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
