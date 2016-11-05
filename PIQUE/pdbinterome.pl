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


my $getfromfile = 0  ;

use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($cB,$fastafile,$idx,$contactfile,$workdir,$p1,$p2,$outfile,$cutoff,$cA,$listfile,$protein);
my ($ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 0 ;
$cutoff = 1 ;
GetOptions(
            "cA=s"=>\$cA ,
            "cB=s"=>\$cB ,
            "getfromfile=i"=>\$getfromfile ,
            "protein=s"=>\$protein ,
            "fastafile=s"=>\$fastafile ,
            "workdir=s"=>\$workdir ,
            "contactfile=s"=>\$contactfile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "idx=i"=>\$idx ,
            "verbose=i"=>\$verbose ,
            "cutoff=i"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);
usage( "Need to give a input file name => option -contactfile ") if(!defined $contactfile);
usage( "Need to give a input file name => option -cA ") if(!defined $cA);
usage( "Need to give a input file name => option -cB ") if(!defined $cB);
usage( "Need to give a input file name => option -workdir ") if(!defined $workdir);

usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);

my @list= util_read_list_sentences($listfile);
my $list = {};
map { s/\s*//g ; $list->{$_} = 1 ; } @list ;


my $infoResidues = {};

### The interacting residues may be provided in a file
my $numberocontacts = 0 ;
if($getfromfile){
	$getfromfile = "$cA.$cB";
	print "Getting interacting residues from $getfromfile\n";
    my $ifh = util_read($getfromfile);
	while(<$ifh>){
        next if(/^\s*$/);
	    next if(/^\s*#/);
		my (@l) = split ;

		foreach my $r (@l){
		    $infoResidues->{$r} = 1 ;
		}
	}
}
else{
  my $ifh = util_read($contactfile);
  while(<$ifh>){
     next if(/^\s*$/);
	 next if(/^\s*#/);

	 if(/EDGE/){
	    my ($x,$A,$B,$H1,$H2,$res1,$res2) = util_ParseContactFile($_);
		if(($A =~ /$cA/ && $B =~ /$cB/)|| ($B =~ /$cA/ && $A =~ /$cB/)){
			# if bothe the chains are the same, add residues from both side
			$numberocontacts++;
		    if($cA eq $cB){
		       $res1 = $res1 . $res2 ;
		    }
			my @res1 = split ",", $res1 ;
			foreach my $r (@res1){
				$infoResidues->{$r} = 1 ;
			}
		     
		}
	 }
   }
   close($ifh);
}



my $resstr = join " ", util_SortResiduesBasedonNumber(3,keys %{$infoResidues});
print "RESIDUES making A=$cA to B=$cB are (getfromfile=$getfromfile,  numberocontacts=$numberocontacts):\n $resstr\n\n";

my $sortprint = {};

## reverse namemap - to get name
my $namemap = {};
if(-e  "namemap"){
	my $ifhnamemap = util_read("namemap");
	while(<$ifhnamemap>){
		next if(/^\s*#/);
		next if(/^\s*$/);
		my ($a,$b) = split ;
		$namemap->{$b} = $a ;
	}
}

## Process each contact file
my $CHAPERONES = {};
foreach my $x (@list){
	my $i = "WORKDIR/$x.contactchains";
	next if($i eq $contactfile);
	my $ifh = util_read($i);
	while(<$ifh>){
	   if(/EDGE/){
	       my ($junk,$A,$B,$H1,$H2,$RES1,$RES2) = util_ParseContactFile($_);
		   my $res1 ;
		   my $theotherone ;
		   if(($A =~ /$cA/ && !($B =~ /$cB/)) || ($B =~ /$cA/ && !($A =~ /$cB/)) ){
		        if($A =~ /$cA/ && !($B =~ /$cB/)){
					$res1 = $RES1;
					$theotherone = $B ;
				}
				else{
					$res1 = $RES2;
					$theotherone = $A ;
				}
			    my @res1 = split ",", $res1 ;

				my $cnt = 0 ;
				my $STR = "";

			    foreach my $r (@res1){
					if(exists $infoResidues->{$r}){
						$cnt++;
						$STR = $STR . " $r ";
					}
				}
				if($cnt > $cutoff){
					my $FILE = "WORKDIR/list.$x";
					die "no list file WORKDIR/list.$x" if(! -e $FILE);
					my ($PDBFILE) = ` cat $FILE  `;
					chomp $PDBFILE ;
					
					my $str = "$cnt ($STR) $x PDB=$PDBFILE   $_ ";
					$theotherone=~ s/_.*//;
					if(!exists $CHAPERONES->{$theotherone} || $cnt > $CHAPERONES->{$theotherone}){
						$CHAPERONES->{$theotherone} = $cnt;
					}
					print " $A $B theotherone =$theotherone \n" if($verbose);
					$sortprint->{$cnt} = [] if(!defined $sortprint->{$cnt});
					push @{$sortprint->{$cnt}}, $str ;

				}
		   }
	   }

	}
}

foreach my $k (sort {$b <=> $a} keys %{$sortprint}){
	my $v = $sortprint->{$k};
	print @{$v} if($verbose);
}

print "### sorted on number of contacts the chaperone makes\n";
foreach my $k (sort {$CHAPERONES->{$b} <=> $CHAPERONES->{$a}} keys %{$CHAPERONES}){
	my $v = $CHAPERONES->{$k};
	print "$k  $v\n";
}



sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
