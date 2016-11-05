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
usage( "Need to give a input file name => option -fastafile ") if(!defined $fastafile);
my ($tablefasta,$infoSeq2PDB,$mapChainedName2Name) = util_parsePDBSEQRESNEW($fastafile,0,0);


my $changed = {};
my @LISTS ;
my @CHLISTS ;
my $N ;
foreach my $k (sort keys %{$tablefasta}){
	my $seq = $tablefasta->{$k};
	
	my $changedseq = ReplaceStringWithEquivalent($seq);
	$changed->{$k} = $changedseq ;

	my @l = split "", $seq ;
	my @cl = split "", $changedseq ;
	if(defined $N){
		die if(@l ne $N);
		die "$changedseq"  if(@cl ne $N);
	}
	else{
		$N = @l ;
	}

	push @LISTS,\@l;
	push @CHLISTS,\@cl;

}

$N = $N -1;

foreach my $idx (0..480){
	my @l ;
    foreach my $L (@CHLISTS){
		my @L = @{$L};
	    my $JJJJ = @L ;
		my $v = $L[$idx];
		push @l, $v ;
	}
	my @uniqL = util_sortuniqArray(@l);
	my $NN = @uniqL ;
	if($NN eq 1){
	    print "$idx(@uniqL) ";
	}
}
print "\n";


sub ReplaceStringWithEquivalent{
	my ($seq) = @_ ;
	$seq =~ s/E/D/g;
	$seq =~ s/K/R/g;
	$seq =~ s/H/R/g;

	$seq =~ s/T/S/g;

	$seq =~ s/Y/F/g;
	$seq =~ s/W/F/g;
	return $seq ;
}

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
