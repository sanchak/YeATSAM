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
use MyConfigs;


use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($duplicate,$blastout,$infile,$p1,$p2,$outfile,$tag,$cutoff,$listfile,$protein);
my ($chooseone,$ignorefile,$donone,@expressions);
my $howmany = 100000 ;
my $verbose = 0 ;
GetOptions(
            "tag=s"=>\$tag ,
            "blastout=s"=>\$blastout ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            #"outfile=s"=>\$outfile ,
            "donone"=>\$donone ,
            "duplicate"=>\$duplicate ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
#usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
#my $ofh = util_write($outfile);

my $orfdir = "FASTADIR_ORF";
usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
usage( "Need to give a tag -option -tag  ") if(!defined $tag);


my @list= util_read_list_sentences($listfile);

my $fhNONE = util_write("$tag.none");
my $fhUNIQ = util_write("$tag.uniq");

PrintALL("ANN/LISTS/trs.list","$tag.list",@list);
PrintALL("ANN/ANNOTATE/trs.a","$tag.anno",@list);
PrintALL("ANN/COMMANDS/trs.c","$tag.commands.csh",@list);
system(" extractindexfromfile.pl -in $tag.anno");

   #my $ifhwarn = util_read("ANN/WARN/$trs.w");
   #my $ifhcommands = util_read("ANN/COMMANDS/$trs.c");
#
   #my $ifhannotate = util_read("ANN/ANNOTATE/$trs.a");


sub PrintALL{
	my ($infilestring,$outfile,@LLL) = @_ ;
	print "Info: Processing $infilestring for $outfile\n";
    my $done = {};
    my $OFH = util_write($outfile);
    foreach my $i (@LLL){
		$i =~ s/ //g;
		my $x = $infilestring ;
		$x =~ s/trs/$i/;
        my $ifhlist = util_read($x);
        while(<$ifhlist>){
   	         next if(/^\s*$/);
		     my ($i) = split ;
			 if(/#none/){
   	              print $fhNONE "$i\n";
				  next ;
			 }
			 if(/#uniq/){
   	              print $fhUNIQ "$i\n";
			 }
		     #warn "$i done $infilestring $outfile in $x "  if(exists $done->{$i} && $i ne "\\cp" && $i ne "fixFastaNames.pl");
		     $done->{$i} = 1 ;
   	         print $OFH $_ ;
        }
        close($ifhlist);
	}
}

system ("wc -l $tag.*");

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
