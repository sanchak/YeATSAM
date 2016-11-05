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
my ($duplicate,$blastout,$infile,$p1,$p2,$outfile,$trs,$cutoff,$unique,$mapfile,$listfile,$protein);
my ($chooseone,$ignorefile,$inorderoflist,$tag,$donone,@expressions);
my $howmany = 100000 ;
my $verbose = 0 ;
my $removeslash = 1 ;
GetOptions(
            "unique=s"=>\$unique ,
            "mapfile=s"=>\$mapfile ,
            "trs=s"=>\$trs ,
            "tag=s"=>\$tag ,
            "blastout=s"=>\$blastout ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "inorderoflist"=>\$inorderoflist ,
            "donone"=>\$donone ,
            "duplicate"=>\$duplicate ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
            "removeslash=i"=>\$removeslash ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);

usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
usage( "Need to give a mapfile -option -mapfile  ") if(!defined $mapfile);
usage( "Need to give a infile -option -infile  ") if(!defined $infile);
$tag = "ttt" if(!defined $tag);
$outfile = "$tag.$infile" if(!defined $outfile);

$outfile =~ s/\///g if($removeslash);


my $ofh = util_write($outfile);
my $ifh = util_read($infile);

my ($table,$XX) = util_maketablefromfile($mapfile);
my $revtable = {};
foreach my $k (keys %{$table}){
	my $v = $table->{$k};
	$revtable->{$v} = $k;
}

my ($table2print,$XXX) = util_maketablefromfile($listfile);

my $ign = 0;
while(<$ifh>){
   my ($trs,@l) = split ;
   if(! exists $revtable->{$trs}){
   	 $ign++;
   	 next ;
   }
   else{
   	  my $v = $revtable->{$trs};
	  if(exists $table2print->{$v}){
	     print $ofh "$v $trs @l\n";
	  }
   }
}
print "$ign = ignored, \n";

system(" wc -l $listfile");
system(" wc -l $outfile");
   
sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
