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
my ($het2pdb,$pdbseqres,$infile,$outfile,$which_tech,$comparewith,$protein);
my (@expressions);
my $size ;
my $verbose = 0 ;
GetOptions(
            "het2pdb=s"=>\$het2pdb ,
            "pdbseqres=s"=>\$pdbseqres ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "comparewith=s"=>\$comparewith ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "size=i"=>\$size ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);

usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a input file name => option -comparewith ") if(!defined $comparewith);
usage( "Need to give a input file name => option -size ") if(!defined $size);


print "READING $infile\n";
my ($infoCompare,$infoSeq2PDBCompare,$mapChainedName2NameCompare) = util_parsePDBSEQRESNEW($comparewith,0);
my ($info,$infoSeq2PDB,$mapChainedName2Name) = util_parsePDBSEQRESNEW($infile,0);

my $done = {};
foreach my $k (keys %{$info}){
	my $seq = $info->{$k} ;
	my $len = length($seq);
    next if($len < $size);

	if($seq =~ /^M/){
		my @l = ($seq =~ /(.)/g);
		my $str = "";
		foreach my $i (1..$size){
			my $idx = $i - 1;
			my $x = $l[$idx];
			$str = $str . $x ;
		}

		if(exists $infoSeq2PDBCompare->{$str}){
			my ($v) = @{$infoSeq2PDBCompare->{$str}};
			print $ofh "$k $str $v\n";
		}
	}

}

system("wc -l $outfile");


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
