#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use MyConfigs;
use MyGeom;
use PDB;
use ConfigPDB;

use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($percentfile,$infile,$outfile,$template,$listfile,$protein);
my (@expressions,$config);
my $aminoacid = 1;
my $force = 0;
my $verbose = 1 ;
my $printTex = 0 ;
my $aadiff = "aadiff.rep";
GetOptions(
            "template=s"=>\$template ,
            "percentfile=s"=>\$percentfile ,
            "printTex=s"=>\$printTex ,
            "aadiff=s"=>\$aadiff ,
            "protein=s"=>\$protein ,
            "config=s"=>\$config ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "aminoacid=i"=>\$aminoacid ,
            "force=i"=>\$force ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_open_or_append($outfile);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a protein pdb id -option -protein  ") if(!defined $protein);

my $CNT = 0 ; 
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my $PWD = cwd;

my ($donetable) = util_maketablefromfile($outfile);
if(exists $donetable->{$protein} && !$force){
	print "Info: protein already processed in $outfile, use force=1\n";
}




my ($string,$firstline) = util_readfasta($infile);

my $mapAA = {}; 

my @ALL = ($string =~ /./g);
my $LEN = @ALL;




my  ($tableThree2One,$tableOne2Three,@sortedSingle) = Config_AACodes();
my $sum = 0; 
print $ofh "$protein $LEN ";

my @lll ;

my $table = {};
$table->{"P"} = 1 ;
$table->{"K"} = 1 ;
$table->{"V"} = 1 ;

my $cntamino = 0 ;
my $sorted = {};
foreach my $one (@sortedSingle){
  	my $i = $one ;
    my @M = ($string =~ /$one/g);
    my $N = @M ; 
	if(exists $table->{$one}){
	   $cntamino = $cntamino+$N;	
	}
	$sum = $sum + $N ; 
	my $per = util_format_float(100*($N/$LEN),1);
	$per =~ s/\s*//g;
    if($aminoacid){
	    #print $ofh " $per ";
	    print $ofh " $one $per ";
		$mapAA->{$one} = $per ; 
	}
	else{
	    print $ofh " $one= $per ";
	}
	push @lll, $per ;
  } 

  if($aminoacid){
     print $ofh "\n";
  }
  else{
     my ($mean,$sd) = util_GetMeanSD(\@lll);
	 my $diff = util_format_float(($lll[0] + $lll[2]) - ($lll[1] + $lll[3]),1) ;
     print  $ofh "sd= $sd diff= $diff \n";
  }

  if($aminoacid && !$cntamino){
  	    warn "Something wrong - maybe you are using nt for amino acids\n";
  }

my $TEMPLATE = {};
if(defined $template){
	($TEMPLATE) = util_maketablefromfile($template);
}


my $ofhdiff = util_open_or_append($aadiff);
  ## print for TEX 
  if($aminoacid){
  	my $TMPFH ;
	if($printTex){
  	    $TMPFH = util_write("$protein.tmp.tex");
	    print $TMPFH "\\addplot [draw=none, fill=green!90] coordinates{ \n";
	}
	my $diff = 0 ;
  	foreach my $k (sort keys %{$mapAA}){
		my $v = $mapAA->{$k} ;
		if(defined $template){
		     my $VTMEP = $TEMPLATE->{$k} ;
			 my $indiviDill = abs($VTMEP - $v);
			 print "$k  abs($VTMEP - $v) = $indiviDill lllllll\n";
			 $diff = $diff + abs(util_format_float($VTMEP - $v,1)) ;
		}
		print $TMPFH "($k, $v)\n" if($printTex);
	}
	print $TMPFH "};\n" if($printTex);
	print $ofhdiff "$protein $diff \n";
  }

  #my $percent = (100 * $sum )/$LEN ; 
  #print $ofh "$protein $percent\n";

  # do stuff with $string

sub usage{
    my ($msg) = @_ ;
	    print $msg , "\n" ;
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
		    die ;
			}
