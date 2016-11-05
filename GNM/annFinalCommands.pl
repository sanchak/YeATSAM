#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use Carp ;
use MyUtils;
use MyGeom;
use PDB;
use ConfigPDB;
use MyGNM;
use Math::Geometry ;
use Math::Geometry::Planar;

## disable perl's warning mechanism
no warnings 'recursion';

use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($mapfile,$isNT,$findcharstring,$forWGS,$checkforsame,$infile,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions,$trs,$blastdir);
my $strict = 1 ;
my $verbose = 0 ;

my $blastcutoff ;
my $percentlength ;
my $percentmatched ;
my $percentidentity ;
my $expectlimit ;

my $postfix = "blast.nt";
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "mapfile=s"=>\$mapfile ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "postfix=s"=>\$postfix ,
            "trs=s"=>\$trs ,
            "checkforsame"=>\$checkforsame ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "blastdir=s"=>\$blastdir ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "blastcutoff=i"=>\$blastcutoff ,
            "strict=i"=>\$strict ,
            "forWGS=i"=>\$forWGS ,
            "verbose=i"=>\$verbose ,
            "percentlength=i"=>\$percentlength ,
            "percentmatched=i"=>\$percentmatched ,
            "percentidentity=i"=>\$percentidentity ,
            "findcharstring=i"=>\$findcharstring ,
            "isNT=i"=>\$isNT ,
            "expectlimit=f"=>\$expectlimit ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a input file name => option -blastcutoff ") if(!defined $blastcutoff);
usage( "Need to give a input file name => option -forWGS ") if(!defined $forWGS);
usage( "Need to give a input file name => option -listfile ") if(!defined $listfile);

my $fname = "LISTS/$outfile/$outfile.$blastcutoff" ;

system ("simpleAppendMatches.pl -inf $fname");
system ("hackAnnoRemoveValues.pl -inf $fname.appended");
if(! $forWGS){
     die "Need map file for lengths when you are grouping" if(!defined $mapfile);
     print "groupBasedonCutoff.pl -outf $fname.group -inf $fname -cutoff $blastcutoff -dir 1\n";
     system ("groupBasedonCutoff.pl -outf $fname.group -inf $fname -cutoff $blastcutoff -dir 1");

     system ("groupSortBasedOnLength.pl -inf $fname.group -map $mapfile");
     
     
     system ("\\rm -f  $fname.GROUPED");
     system("twolists.pl $fname.anno $fname.group.allingroup -ver 0 ");
     system("cat ofhinAbutnotinB  $fname.group.first.len.sort.first> $fname.GROUPED");
}




system ("mappingAddCount.pl -in $fname -outf $fname.found -ignoresingle");
system ("\\rm -f  $fname.anno.real");
system ("cat $fname.anno | grep -v Warning > $fname.anno.real");
system ("extractindexfromfile.pl -in $fname.anno.real");
system ("extractindexfromfile.pl -in $fname.anno.real -out $fname.list");
system ("wc -l $fname.*");
system("twolists.pl $listfile $fname.anno.real -ver 0 ");
system ("\\mv -f ofhinAbutnotinB $fname.anno.no");


system ("mappingAddCount.pl -infile $fname.appended.removevalues -ignoresingle");
system ("sortOnLast.pl -infile $fname.anno.real -idx 1 -rever ");
system ("specialCreateTexTableForAnno.pl -infile $fname.anno.real.sort");

print "IN CASE YOU WANT TO GROUP: groupBasedonCutoff.pl -in $fname.appended.removevalues.mapped.pw -out jjj -cutoff 0\n";




sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
