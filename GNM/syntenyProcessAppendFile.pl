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
my ($fastafile,$idx,$postfix,$p1,$p2,$cutoff,$which_tech,$listfile,$protein);
my ($significance,$ignorefile,@expressions);
my $howmany = 100000 ;
$cutoff = 5 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "fastafile=s"=>\$fastafile ,
            "postfix=s"=>\$postfix ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "idx=i"=>\$idx ,
            "significance=i"=>\$significance ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a input file name => option -listfile ") if(!defined $listfile);
usage( "Need to give a input file name => option -significance ") if(!defined $significance);

my @list= util_read_list_sentences($listfile);
my $list = {};
map { s/\s*//g ; $list->{$_} = 1 ; } @list ;

my $ONLYONEDIR = "ONLYONE.$significance";
my $SYNTENY = "SYNTENY.$significance";
my $MEANSD = "MEANSD";

system("mkdir -p $MEANSD");
system("rm -f $MEANSD/*");

system("mkdir -p $ONLYONEDIR");
system("rm -f $ONLYONEDIR/*");

system("mkdir -p $SYNTENY");
system("rm -f $SYNTENY/*");

print "Info: Using cutoff $cutoff\n";

my @infos ;
foreach my $i (@list){
   my $infile = "PARSEBLASTOUT/MMM." . $i .  ".$significance.appended";
   my ($info,$sorted) = ProcessOneAppendFile($infile);

   my $ofhsynteny = util_open_or_append("$SYNTENY/$i.out");
   if(1){
   	   my $ToSort ={};
       foreach my $k ( sort { $sorted->{$b} <=> $sorted->{$a}} keys %{$info}){
	       my @l = @{$info->{$k}};

		   ## no point of synteny if match is only 1
		   next if(@l < 2);
		   my $first = $l[0];
   	       my $str = "$k @l";
		   $ToSort->{$str} = $first ;
       }

	   foreach my $k (sort {$ToSort->{$a} <=> $ToSort->{$b}} keys %{$ToSort}){
	   	   my @lll = split " ", $k ;
		   my $name = shift @lll ;
		   my $NUM = @lll ;
		   my $s = shift @lll ;
		   my $e = $lll[@lll -1];
		   my $diff = abs($e -$s);
		   #my $average = util_format_float($diff/($NUM*$NUM),1); ## NUM sq - give some value to more
		   my $average = util_format_float($diff/($NUM),1); 

		   printf $ofhsynteny "%40s %6s %5s %4s %3s %5s %5s\n", $name,$i,$s,$e,$NUM,$diff,$average ;
	   }
   }
   close($ofhsynteny);
   push @infos, $info ;
}

my $N = @infos - 1  ;

my $ofhdiff = util_write("diffNums");

my $ALLSCAFFSINCONSIDERATION = {};
foreach my $i (0..$N){


	## counts the number of occurences in others 
	my $COUNTOccurancesinOthers = {};

    foreach my $j (0..$N){
		next if($i eq $j);
	    my $other = $infos[$j];
		foreach my $k (keys %{$other}){
			$ALLSCAFFSINCONSIDERATION->{$k} = 1 ;
			my $NUM = @{$other->{$k}};
			$COUNTOccurancesinOthers->{$k} = 0 if(! defined $COUNTOccurancesinOthers->{$k});
			$COUNTOccurancesinOthers->{$k} = $COUNTOccurancesinOthers->{$k} + $NUM ;
		}
    }

	my $infocus = $infos[$i];
	my $chrname = $list[$i];
	foreach my $k (sort keys %{$infocus}){
		$ALLSCAFFSINCONSIDERATION->{$k} = 1 ;
	    my $NUM = @{$infocus->{$k}};
		if(! exists $COUNTOccurancesinOthers->{$k}){
            my $ofhonlyone = util_open_or_append("$ONLYONEDIR/$chrname");
			print $ofhonlyone "$k $NUM \n";
		}
		else{
			my $NOTHERS = $COUNTOccurancesinOthers->{$k} ;
			my $difff = $NUM - $NOTHERS ;
			print $ofhdiff "$k with $difff difff in $chrname....\n";
		}
	}
}


## This counts the number of occurences of each scaffold in each Chromosome
## by sorting, you can actually get to see which are similar
my $ofhNICEMAP = util_write("nicemap.$significance.csv");
my $ofhMeanSD = util_write("$MEANSD/meansd.$significance.csv");
my $NSCAFF = (keys %{$ALLSCAFFSINCONSIDERATION});
print "There are $NSCAFF mapped scaffolds\n";
foreach my $k (sort keys %{$ALLSCAFFSINCONSIDERATION}){
	print $ofhNICEMAP "$k,  ";

	my @values ; 
    foreach my $i (0..$N){
	    my $infocus = $infos[$i];
	    my $NUM = 0 ;
		if(exists $infocus->{$k}){
		   $NUM = @{$infocus->{$k}};
		}	
		push @values, $NUM ;
		print $ofhNICEMAP " $NUM , ";
	}
	print $ofhNICEMAP " \n";
    my ($mean,$sd) = util_GetMeanSD(\@values);

	#my $meanmultsd = $mean * $sd ;
	my $meanmultsd = $mean / $sd ;
	printf $ofhMeanSD "%40s ", $k ;
	print $ofhMeanSD " @values $mean $sd $meanmultsd\n";



}
close($ofhNICEMAP);


my $ofhtmp = util_write("tttt");
foreach my $i (@list){
	print $ofhtmp  "sortOnLast.pl -rev -in $ONLYONEDIR/$i  & \n";
	print $ofhtmp  "sortOnLast.pl  -in $SYNTENY/$i.out  & \n";
	print $ofhtmp  "sortOnLast.pl  -in $MEANSD/meansd.$significance.csv  & \n";
}
close($ofhtmp);
system("runInBack.csh tttt");

sleep(1);
system("wc -l $ONLYONEDIR/chr*sort");
unlink "tttt";

sub ProcessOneAppendFile{
  my ($infile) = @_ ;
  my $outfile =  "$infile.synout";
  my $ofh = util_write($outfile);
  my $ifh = util_read($infile);
  my $cnt = 0 ;
  my $info = {};
  while(<$ifh>){
       next if(/^\s*$/);
	 next if(/^\s*#/);

	 ## note this is order dependent...
	 $cnt++;

     next if(/0 # only one/);
	 my ($nm,@l) = split ; 
	 while(@l){
	 	my $scaff = shift @l ;
	 	shift @l ;
		next if($scaff =~ /MITO/ || $scaff =~ /CHLORO/);
		if(! exists $info->{$scaff}){
			 $info->{$scaff} = [];
		}
		push @{$info->{$scaff}} , $cnt ;
	 }
  }
  close($ifh);

  my $sorted = {};
  foreach my $k (keys %{$info}){
	 my @l = @{$info->{$k}};
	 my $N = @l ;
     $sorted->{$k} = $N ;
  }

  my $inforet = {};
  foreach my $k (sort { $sorted->{$b} <=> $sorted->{$a} } keys %{$info}){
      my $v = $sorted->{$k};
	  if($v > $cutoff){
	  	  $inforet->{$k} = $info->{$k} ;
  	      print $ofh "$k $v\n";
	  }
  }


  close($ofh);
  return ($inforet,$sorted);
}



sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
