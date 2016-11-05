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
my ($fastadir,$map2trs,$infile,$p1,$p2,$orfdir,$outfile,$cutoff,$listfile,$orf,$trs);
my ($fastafile,@expressions);
my $howmany = 100000 ;
my $verbose = 0 ;
GetOptions(
            "listfile=s"=>\$listfile ,
            "trs=s"=>\$trs ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "map2trs=s"=>\$map2trs ,
            "p2=s"=>\$p2 ,
            "orf=s"=>\$orf ,
            "fastadir=s"=>\$fastadir ,
            "outfile=s"=>\$outfile ,
            "orfdir=s"=>\$orfdir ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a output file name => option -orfdir ") if(!defined $orfdir);
my $ofh = util_write($outfile);
usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
usage( "Need to give a fastadir -option -fastadir  ") if(!defined $fastadir);
usage( "Need to give a map2trs -option -map2trs  ") if(!defined $map2trs);

my $info = {};

my $IFH = util_read($listfile);
while(<$IFH>){
	next if(/^\s*$/);
	next if(/^\s*#/);
	my (@l) = split ; 
	my $name = $l[0];
    my $fhlist = util_write("TMPDIR/$name.list");
	$, = "\n";
	print $fhlist "@l \n";
	close($fhlist);

    my $tmplist = "$name.ttttt";

	print $ofh "unlink $tmplist \n";
	print $ofh "touch $tmplist\n";
	foreach my $l (@l){
		print $ofh "grep -w $l $map2trs >> $tmplist\n";
	}
	print $ofh "extractCodonFromORFs.pl -outf codonbias -lis $tmplist -cutoff 40 -fasta $fastadir -orfdir $orfdir > ! $name.log  \n";
}

print "Writing commands to $outfile\n";

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
