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
my ($duplicate,$blastout,$infile,$p1,$p2,$outfile,$trs,$cutoff,$unique,$listfile,$protein);
my ($chooseone,$ignorefile,$inorderoflist,$tag,$donone,@expressions);
my $howmany = 100000 ;
my $ignorePostfix ;
my $verbose = 0 ;
my $removeslash = 1 ;
GetOptions(
            "unique=s"=>\$unique ,
            "trs=s"=>\$trs ,
            "ignorePostfix"=>\$ignorePostfix ,
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
usage( "Need to give a infile -option -infile  ") if(!defined $infile);
$tag = "ttt" if(!defined $tag);
$removeslash  = 0 if (defined $outfile);
$outfile = "$tag.$infile" if(!defined $outfile);

$outfile =~ s/\///g if($removeslash);


my $ofh = util_write($outfile);
my $ofhmissed = util_write("missed");
my $ifh = util_read($infile);

my $table = {};
my @LLL= util_read_list_words($listfile);
map {  my @ll = split ;  
my $nm = $ll[0];

#$nm =~ s/\.MER.*//;

if(defined $ignorePostfix){
$nm =~ s/_A//;
$nm =~ s/_B//;
$nm =~ s/_C//;
}

$table->{$nm} = 1 ; } @LLL ;
my @list = (keys %{$table});

my $cnt = 0 ;
my $missed = 0 ;
my $printed = {};
while(<$ifh>){

   my $orog = $_;

   # just hardcoded, essentially we take the first
   s/ln -s//;
   #s/.ORF.*//;
   my ($trs) = split ;
   if(exists $table->{$trs}){
   	   
	   next if(defined $unique && exists $printed->{$trs});
   	   $cnt++;
   	   print $ofh "$orog" if(!defined $inorderoflist);
	   $printed->{$trs} = $orog;
   }
}

if(defined $inorderoflist){
	foreach my $i (@LLL){
	   if(exists $printed->{$i}){
	      my $orog = $printed->{$i};
   	      print $ofh "$orog" ;
	   }
	}
}

foreach my $k (keys %{$table}){
	if(! exists $printed->{$k}){
		#$k =~ s/_A//;
		#$k =~ s/_B//;
		#$k =~ s/_C//;
	}
	if(! exists $printed->{$k}){
		$missed++;
		print $ofhmissed "$k\n";
	}
}

print "Wrote $cnt to $outfile , missed $missed\n";
   
sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
