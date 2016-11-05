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
my ($blastdir,$fdirDB,$fdirTRS,$idx,$infile,$postfix,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "blastdir=s"=>\$blastdir ,
            "fdirDB=s"=>\$fdirDB ,
            "fdirTRS=s"=>\$fdirTRS ,
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
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);
usage( "Need to give a input file name => option -fdirDB ") if(!defined $fdirDB);
usage( "Need to give a input file name => option -fdirTRS ") if(!defined $fdirTRS);
usage( "Need to give a input file name => option -blastdir ") if(!defined $blastdir);


my $info = {};

## collate - as there can be multiple
#  This reads ".20.mapone" - which stores the match with the numbers
#  in the non-masked version, there will be only one entry

my $MASKEDVERSION = 0 ;
while(<$ifh>){
     next if(/^\s*$/);
	 next if(/^\s*#/);
	 my ($nm,@l) = split ; 
	 if(!defined $info->{$nm}){
		$info->{$nm} = [];
	}
	else{
		$MASKEDVERSION = 1 ;
	}
	push @{$info->{$nm}}, @l ;
}
close($ifh);


## Find the match with the highest numbers...
my $newinfo = {};
foreach my $k (keys %{$info}){
	my @l = @{$info->{$k}};
	my $bestmatch ;
	my $bestnumber ;
	if($MASKEDVERSION){
	 while(@l){
		my $match = shift @l ;
		my $number = shift @l ;
		if(!defined $bestnumber){
			$bestmatch = $match ;
			$bestnumber = $number ;
		}
		else{
			if($bestmatch eq $match){
				$bestnumber = $bestnumber + $number ;
			}
			else{
				if($match > $bestmatch){
			       $bestmatch = $match ;
			       $bestnumber = $number ;
				}
			}
		}
	 }
	}
	else{
			die if(@l ne 2);
			$bestmatch = shift @l ;
			$bestnumber = shift @l ;
	}
	$newinfo->{$k} = $bestmatch ;
	#print "$k has bestmatch $bestmatch with bestnumber = $bestnumber \n";
}

my $ofh = util_write($outfile);
system("mkdir -p $blastdir");
foreach my $k (keys %{$newinfo}){
	 	my $i = $newinfo->{$k} ;
		## pairwise ...
	 	print $ofh "blastp -query $fdirTRS/$k.ALL.1.fasta -subject $fdirDB/$i.ALL.1.fasta -out $blastdir/$k.blast.nt \n"; 
}


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
