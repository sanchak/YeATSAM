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

usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
my $CNT = 0 ; 
 my ($RESULTDIR,$PDBDIR) = util_SetEnvVars();

# for reading fasta
#my ($fe,$x) =  util_readfasta($infile);
#chomp $x ;
my $PWD = cwd;


my $ignoretable ;
if(defined $ignorefile){
    ($ignoretable) = util_maketablefromfile($ignorefile);
}

my  ($seconds, $microseconds) = gettimeofday;

my @list= util_read_list_sentences($listfile);
my $info = {};
foreach my $protein (@list){
    my $pdb = "$PDBDIR/$protein.pdb";

   my $ifh = util_read($pdb);
   while(<$ifh>){
     next if(/^\s*$/);
	 next if(/^\s*#/);

     chomp ;
	 if(/ORGANISM_SCIENTIFIC: /){
	 	s/.*ORGANISM_SCIENTIFIC://;
	 	s/;//;
	 	s/ //g;
		my $organism = $_ ;
		$info->{$organism} = [] if(!defined $info->{$organism});
		push @{$info->{$organism}}, $protein ;
	 }
   }
}

foreach my $organism (sort keys %{$info}){
	my @pdbs = @{$info->{$organism}};
	my $N = @pdbs ;
	print $ofh "$organism $N  @pdbs \n";
}

#my ($N,$first,$last,$range,$mean,$sd) = util_GetStatOfList(@l);
#my ($table) = util_make_list(\@list);
#my ($tablemerge,$newN) = util_table_merge($t1,$t2);
#my ($table) = util_mapFullLinetofirst($fname);
#my ($table,$N) = util_maketablefromfile($fname);
#my ($common,$inAbutnotinB,$inBbutnotinA) = util_table_diff($t1,$t2);
#my $junk = util_writelist2file($fname,@list);


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
