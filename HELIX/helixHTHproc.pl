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
my ($cutofflengthofturn,$infile,$outfile,$printPDB,$which_tech,$listfile,$protein);
my (@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
$cutofflengthofturn = 4 ;
print STDERR "Info: cutofflengthofturn = $cutofflengthofturn\n";
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "printPDB"=>\$printPDB,
            "howmany=i"=>\$howmany ,
            "verbose=i"=>\$verbose ,
            "cutofflengthofturn=i"=>\$cutofflengthofturn ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);
usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
usage( "Need to give a protein pdb id -option -protein  ") if(!defined $protein);
my $CNT = 0 ; 
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();
my $PWD = cwd;

my  ($seconds, $microseconds) = gettimeofday;


my $info = {};
while(<$ifh>){
	chop ;
	 my ($pdb,$helix) = split ; 
	 my $key = "$pdb.$helix";
	 $info->{$key} = $_ ; 
}
close($ifh);


my @LLL= util_read_list_sentences($listfile);

## helix in order  +ve and then -ve
Run(1,@LLL);
## helix in reverse order  -ve and then +ve
Run(0,@LLL);




sub Run{

   my ($what,@list) = @_ ;
   foreach my $l (@list){
	  my ($h1,$h2,$prot,$lengthofturn,$start,$end) = split " ", $l;
	  next if(!($h1 =~ /HELIX/));

      if(defined $printPDB){
	      my $tableofres = {};
		  foreach my $i ($start..$end){
		      $tableofres->{$i} = 1 ;
		      }

           my $pdb = "$PDBDIR/$prot.pdb";
           my $pdb1 = new PDB();
           $pdb1->ReadPDB($pdb);
		   	my $OOO = "$h1.$h2.pdb";
		    my $ofhhelix= util_write($OOO);
			print "Writing file $OOO \n";
			$pdb1->ReadPDBAndWriteSubset($pdb,$tableofres,$ofhhelix);
			close($ofhhelix);

	}


	my $donext = 0 ; 

	my $key1 = "$prot.$h1";
	my $key2 = "$prot.$h2";
	if (!exists $info->{$key1}){
		print "$key1 does not exist in allvalues \n";
		$donext = 1 ;
	}
	if (!exists $info->{$key2}){
		print "$key2 does not exist in allvalues \n";
		$donext = 1 ;
	}
	next if($donext);

	#next if($info->{$h1} eq -1 || $info->{$h2} eq -1);

	my $a = $info->{$key1};
	my $b = $info->{$key2};

	my @l1 = split " ",  $a;
	my @l2 = split  " ", $b;

	shift @l1 ; shift @l1 ; shift @l1 ;
	shift @l2 ; shift @l2 ; shift @l2 ;


	if($what eq 1){
	    
		my $IsPosAH = IsPosAH(@l1);
		next if(!$IsPosAH);
		my $IsHydAH = IsNegAH(@l2);
		next if(!$IsHydAH);
	}
	else{
		my $IsHydAH = IsNegAH(@l1);
		next if(!$IsHydAH);
		my $IsPosAH = IsPosAH(@l2);
		next if(!$IsPosAH);
	}

	my $totallen = $lengthofturn + $l1[1] + $l2[1];

	next if ($lengthofturn > $cutofflengthofturn);
	#next if ($totallen > 30);

	print $ofh  "$a \t\t $b turn=$lengthofturn totallen=$totallen\n";

  }
  print $ofh  "\n";
}



sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
sub IsPosAH{
	my (@l1) = @_ ;
	return 0 if($l1[2] < 8);
	return 0 if($l1[3] < 0.8);
	return 0 if($l1[4] < 5);
	print "@l1 is +ve \n" if($verbose);
	return 1 ;
}
sub IsNegAH{
	my (@l1) = @_ ;
	return 0 if($l1[2] < 8);
	return 0 if($l1[3] > 0.3);
	return 0 if($l1[4] < 4);
	print "@l1 is -ve \n" if($verbose);
	return 1 ;
}

sub IsHydAH{
	my (@l2) = @_ ;
	return 0 if($l2[1] < 6 );
	return 0 if($l2[2] > 4);
	return 0 if($l2[4] > 3);
	print "@l2 is hyd \n" if($verbose);
	return 1 ;
}
