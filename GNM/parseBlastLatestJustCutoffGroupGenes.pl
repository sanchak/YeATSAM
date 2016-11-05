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
my ($checkforsame,$infile,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions,$trs);
my $howmany = 100000 ;
my $verbose = 1 ;

my $percentlength = 10;
my $percentmatched = 70;
my $percentidentity = 30;
$cutoff = 0.00000000001 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "trs=s"=>\$trs ,
            "checkforsame"=>\$checkforsame ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "verbose=i"=>\$verbose ,
            "percentlength=i"=>\$percentlength ,
            "percentmatched=i"=>\$percentmatched ,
            "percentidentity=i"=>\$percentidentity ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_open_or_append($outfile);


usage( "Need to give a input file name => option -listfile ") if(!defined $listfile);

my $ifh = util_read($listfile);
my $ID2TRS = {};
while(<$ifh>){

   my ($trs,$infile) = split ;
   my ($info,$querylength) = util_PARSEBLAST($infile);

   foreach my $k (@{$info}){
	   my $org = $k ;
       my ($name,$description,$blastscore,$subjectlength,$iden_percent,$subjectmatched,$querystart,$queryend,$subjectstart,$subjectend,$expect) = split " ",$k ;
	   if(($description =~ /(chromosome|hypothe|unnamed|uncharacterized)/i)){
	     next ;
	   }
	   if($expect < $cutoff){ 
				 my ($realname) = ($k =~ /\|(.*)\|/);
				 print "$realname\n";
				 if(!defined $ID2TRS->{$realname}){
				 	$ID2TRS->{$realname} = [];
				 }
				 push @{$ID2TRS->{$realname}}, $trs ;
	   }
   }
}

my $done = {};
my $UNIQUE = {};
my $OFH = util_write("pwgroups");
foreach my $realname (keys %{$ID2TRS}){
	my @sl = sort @{$ID2TRS->{$realname}} ;


	{
		my @t = @sl ;
		while(@t){
			my $a = shift @t ;
			foreach my $b (@t){
				print $OFH "$a $b 0 \n";
			}
		}

	}
	my $key = join @sl , "" ;
	print "$key \n";
	if(!exists $done->{$key}){
		$done->{$key} = 1 ;
		$UNIQUE->{$realname} = $ID2TRS->{$realname} ;
		print $ofh "$key $realname \n";
	}
}



sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
