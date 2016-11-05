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
my $verbose = 0 ;
$cutofflengthofturn = 4 ;
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
print STDERR "Info: cutofflengthofturn = $cutofflengthofturn\n";

my  ($seconds, $microseconds) = gettimeofday;

my $helices = {};
my $specifichelix = {};
my $processed = 0 ;
my $notprocessed = 0 ;
my $errorfh = util_write("errors");
while(<$ifh>){
     my ($pdb, $helixid, $s,$e,$len,$Lhydro,$PERCENTPOS,$TOTALCHARGED) = split ;
	 my @list = split ;
	 my $N = @list ;
	 if($N ne 8){
	 	$notprocessed++ ;
		print $errorfh "$_";
		next ;
	 }
	 else{
	 	$processed++;
	 }
	 $helices->{$pdb} = [] if(!defined $helices->{$pdb} );
	 push @{$helices->{$pdb}}, \@list ;
	 $specifichelix->{"$pdb.$helixid"} = \@list;
}
print "Info: processed = $processed, notprocessed = $notprocessed \n";



foreach my $protein (sort keys %{$helices}){
	print "Processing $protein\n" if($verbose);
    my $info = {};
    my $prevend ;
    my $helixidprev ;
    my $joiningstretch = -1 ;
    
    
    my $IsPosAHTable = {};
    my $IsNegAHTable = {};
    my $IsHydAHTable = {};
    my $IsConnectedToPrev = {};
    my $IsConnectedToNext = {};
    my @IDS ;
    foreach my $line (@{$helices->{$protein}}){
    	 ## fix this shifting...
    	 my @l1 = @{$line} ;
         my ($pdb, $helixid, $s,$e,$len,$Lhydro,$PERCENTPOS,$TOTALCHARGED) = @l1 ;
    	 shift @l1 ; shift @l1 ; shift @l1 ; 
    
    	 if($helixid =~ /HELIX/ && $pdb eq $protein){
    	      push @IDS, $helixid ;
    	 	  my $IsPosAH = IsPosAH(@l1);
    	 	  my $IsNegAH = IsNegAH(@l1);
    	 	  my $IsHydAH = IsHydAH(@l1);
    		  my $LEN = $e - $s ;
    		  if(defined $prevend){
    		  	   $joiningstretch = $s - $prevend;
    		  }
    		  else{
    		  	   $prevend = $e ;
    		  }
    
    		  $IsPosAHTable->{$helixid} = 1 if($IsPosAH);
    		  $IsNegAHTable->{$helixid} = 1 if($IsNegAH);
    		  $IsHydAHTable->{$helixid} = 1 if($IsHydAH);
    
    		  $info->{$helixid}->{IsPosAH} = $IsPosAH ;
    		  $info->{$helixid}->{IsNegAH} = $IsNegAH ;
    		  $info->{$helixid}->{IsHydAH} = $IsHydAH ;
    		  $info->{$helixid}->{LEN} = $LEN ;
    		  
    		  if($joiningstretch > 0 && $joiningstretch < $cutofflengthofturn){
			  	 #print "$joiningstretch $helixid $helixidprev \n";
    		  	  $IsConnectedToNext->{$helixidprev} = $helixid ;
    		  	  $IsConnectedToPrev->{$helixid} = $helixidprev ;
    		  }
    	      $helixidprev = $helixid ;
    
    		  print "\t$IsPosAH $IsNegAH $IsHydAH $LEN $joiningstretch ...\n" if($verbose);
    	 }
    }
    
    ## Now write the algorithms ...
	## 3Helices - nothing found, hmmm
	if(1){
    foreach my $helixid (@IDS){
    	if(exists $IsConnectedToNext->{$helixid} && exists $IsConnectedToPrev->{$helixid}){
    		my $previd = $IsConnectedToPrev->{$helixid} ;
    		my $nextid = $IsConnectedToNext->{$helixid} ;
    		#if($IsPosAHTable->{$helixid}){
    			#if($IsHydAHTable->{$previd}){
    			    #if($IsHydAHTable->{$nextid}){
    					print $ofh "3HELICES $protein $previd $helixid $nextid \n";
    				#}
    			#}
    		#}
    	}
    }
	}


	if(0){
	my $doneonedir = {};
    foreach my $helixid (@IDS){
    	if(exists $IsConnectedToNext->{$helixid} || exists $IsConnectedToPrev->{$helixid}){
    		my $otherid = exists $IsConnectedToPrev->{$helixid} ? $IsConnectedToPrev->{$helixid}: $IsConnectedToNext->{$helixid} ;
    		if(($IsPosAHTable->{$helixid} && $IsHydAHTable->{$otherid}) 
			       || ($IsPosAHTable->{$otherid} && $IsHydAHTable->{$helixid})){

				my $id1 = "$protein.$helixid";
				my $id2 = "$protein.$otherid";
				if(!exists $doneonedir->{$id1}){
					$doneonedir->{$id1} = 1;
					$doneonedir->{$id2} = 1;
				}
				else{
					next ;
				}
				my $line1 = $specifichelix->{$id1};
				my $line2 = $specifichelix->{$id2};
				my @l1 = @{$line1};
				my @l2 = @{$line2};
    			print $ofh "2HELICES @l1 ..... @l2 ... $protein $helixid $otherid \n";
    		}
    	}
    }
	}

	if(0){
    foreach my $helixid (@IDS){
    		if($IsPosAHTable->{$helixid}){
				    my $line = $specifichelix->{"$protein.$helixid"};
					my @l = @{$line};
					$, = " , " ;
    				print $ofh "SINGLE @l    \n";
			 }
    }
	}
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
	my (@l1) = @_ ;
	return 0 if($l1[1] < 6 );
	return 0 if($l1[2] > 4);
	return 0 if($l1[4] > 3);
	print "@l1 is hyd \n" if($verbose);
	return 1 ;
}
