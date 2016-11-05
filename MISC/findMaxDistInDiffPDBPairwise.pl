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



use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($tag,$fastafile,$idx,$infile,$ahbsfile,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions);
my $cntmax = 10 ;
my $verbose = 1 ;
my $removefile = 0 ;
$cutoff = 3.7 ;
my $pairwise ;
my $findcommon = 0 ;
my $unmappedisDNA = 1;
my $pymolcommands = 1;
print "Info:unmappedisDNA =$unmappedisDNA\n";
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "tag=s"=>\$tag ,
            "protein=s"=>\$protein ,
            "fastafile=s"=>\$fastafile ,
            "ahbsfile=s"=>\$ahbsfile ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "cntmax=i"=>\$cntmax ,
            "removefile=i"=>\$removefile ,
            "idx=i"=>\$idx ,
            "pymolcommands=i"=>\$pymolcommands ,
            "findcommon=i"=>\$findcommon ,
            "pairwise=i"=>\$pairwise ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
#die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -tag ") if(!defined $tag);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a output file name => option -cutoff ") if(!defined $cutoff);
usage( "Need to give a output file name => option -pairwise ") if(!defined $pairwise);
my $ofh = util_write($outfile);

my @list ;
if(defined $listfile){
	print "Reading $listfile\n";
    @list= util_read_list_sentences($listfile);
}
else{
	@list = @ARGV ;
}
 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();


my $OUTDIR = "MAXDIST/$tag";
system ("mkdir -p $OUTDIR");

my $namemap = {};
if(-e  "namemap"){
	my $ifhnamemap = util_read("namemap");
	while(<$ifhnamemap>){
		next if(/^\s*#/);
		next if(/^\s*$/);
		my ($a,$b) = split ;
		$namemap->{$a} = $b ;
	}
}



## You can provide information about helices/betasheets
my $EXTRARG =  "" ;
if(defined $ahbsfile){
		$EXTRARG = " -ah $ahbsfile" ;
}



my $info = {};
my @FORFINDINGCOMMON ;


## pairwise=1 takes a list and compares pairwise - this is done for finding which PDBs are in contact:w
## Otherwise, pairwise=0 we expect a list of pairs - each pair is compared to each other. 
##			This is used to see how different proteins (MCM2, SPT, DAXX) ligand a given protein (H3)...		 
##			The protein in consideration (H3) should be first
##          There are two options here 
##				findcommon=1
##				findcommon=0 - when doing this it is required to look at pairs in both sense

my @SECONDPDBS ;
if($pairwise){
	if(@list < 2){
		print "Nothing to do ";
		exit ;
	}
    my $iter = combinations(\@list, 2);
    while (my $c = $iter->next) {
	    my ($x,$y) = @{$c} ;
	    ProcessPair($x,$y,$EXTRARG);
    }
}
else{
	my $ifh= util_read($listfile);
	$cntmax = 1000 ;
	$EXTRARG = "  -cutoff $cutoff ";
	while(<$ifh>){
		next if(/^\s*#/);
		next if(/^\s*$/);
	    my ($x,$y) = split ;
		push @SECONDPDBS,$y ;
	    ProcessPair($x,$y,$EXTRARG);
	}
}


print "### find common in all .... \n";
my $TABLE1 = shift @FORFINDINGCOMMON ;
foreach my $k (sort {$a <=> $b} keys %{$TABLE1}){

	my $V1 = $TABLE1->{$k};
	my $V2 ;
	my $inall = 1 ;
	foreach my $t (@FORFINDINGCOMMON){
		if(! exists $t->{$k}){
			$inall = 0 ;
			last ;
		}
		else{
	        $V2 = $t->{$k};
		}
	}
	if($inall){
		print " Common in all $k $V1 $V2 \n";
	}
}
print "### find common in all ENDS .... \n";


if($pairwise eq 0){

	my @pairs ;
	## for common, we dont need to reverse a pair - for difference we do 
	#if($findcommon){
	if(0){
	    my $j = shift @SECONDPDBS ;
	    foreach my $k (@SECONDPDBS){
		    die "Cant be same" if($j eq $k);
			push @pairs, $j;
			push @pairs, $k;
		}
	}
	else{
	     foreach my $j (@SECONDPDBS){
	         foreach my $k (@SECONDPDBS){
		         next if($j eq $k);
			     push @pairs, $j;
			     push @pairs, $k;
		     }
		}
	}

	while(@pairs){
		my $j = shift @pairs;
		my $k = shift @pairs;
		my $tabj = $info->{$j} ;
        {
		   next if($j eq $k);
		   my $tabk = $info->{$k} ;
	       my $origa = $j ;
	       my $origb = $k ;
		   if($unmappedisDNA){

	           my $IsDNAA = IsDNA($j);
	           my $IsDNAB = IsDNA($k);
			   die if(exists $namemap->{$j} && $namemap->{$j} ne "DNA" &&  $IsDNAA);
			   die if(exists $namemap->{$k} && $namemap->{$k} ne "DNA" && $IsDNAB);

	           $namemap->{$j} = "DNA" if($IsDNAA && !exists $namemap->{$j});
	           $namemap->{$k} = "DNA" if($IsDNAB && !exists $namemap->{$k});
		   }

	       $j = $namemap->{$j} if(exists $namemap->{$j});
	       $k = $namemap->{$k} if(exists $namemap->{$k});
		   print " ======= PROCESSIN pairwise=$pairwise findcommon=$findcommon $j($origa)  $k($origb) ========= \n";
		   my @MISSINGRESIDUES ;
		   my @MISSINGATOMSINANOTINB ;
		   my @MISSINGATOMSINBNOTINA ;
		   foreach my $key (sort {$a <=> $b}  keys %{$tabj}){
		   	   if(!exists $tabk->{$key}){
			   	   my $Tj = $tabj->{$key} ;
				   my ($justone) = (keys %{$Tj});

				   if(!$findcommon){
			   	       push @MISSINGRESIDUES,$justone ;
				   }
			   }
			   else{
			   	   my $Tj = $tabj->{$key} ;
			   	   my $Tk = $tabk->{$key} ;
				   my ($common,$inAbutnotinB,$inBbutnotinA) =  util_table_diff($Tj,$Tk);
				   my $commonN = @{$common};

				   if($findcommon){
				      if($commonN){
					  	  my ($justone) = @{$common} ;
				   	      print " $justone ";
					  }
					  next ;

				   }
					my @ignore = qw (CE1 CE CG CD1 CD2 C O NE CA CB CZ N CG2 CD);
					my $ignore = util_make_table(\@ignore);
				   	my @lll1 ; 
					foreach my $a1 (@{$inAbutnotinB}){
			               my ($res1,$num1,$type1) = split "/", $a1 ;
					   	   if(! exists $ignore->{$type1}){
						   	push @lll1, $a1 ;
						   }
					}
				   	my @lll2 ; 
					foreach my $a1 (@{$inBbutnotinA}){
			               my ($res1,$num1,$type1) = split "/", $a1 ;
					   	   if(! exists $ignore->{$type1}){
						   	push @lll2, $a1 ;
						   }
					}
				   my $inAbutnotinB_N = @lll1;
				   my $inBbutnotinA_N = @lll2;
				   $, =  "\t" ;


					foreach my $x (@lll1){
					   	   push @MISSINGATOMSINANOTINB,$x ;
					}
					foreach my $x (@lll2){
					   	   push @MISSINGATOMSINBNOTINA,$x ;
					}
			   }
		   }
		   if(!$findcommon){
		   	 my $N1 = @MISSINGRESIDUES;
		   	 my $N2 = @MISSINGATOMSINANOTINB;
		   	 my $N3 = @MISSINGATOMSINBNOTINA;
			 if($N1){
			 	  print "MISSINGRESIDUES @MISSINGRESIDUES\n";
			 }
			 if($N2){
			 	  print "MISSINGATOMSINANOTINB @MISSINGATOMSINANOTINB\n";
			 }
			 if($N3){
			 	  print "MISSINGATOMSINBNOTINA @MISSINGATOMSINBNOTINA\n";
			 }
		   }
	    }
	    print "\n";
    }
}

sub ProcessPair{
	my ($A,$B,$extraarg) = @_ ;
	#print "$A $B ......\n";
	my $infile = "$OUTDIR/$A.$B.maxdist.out";
	if(! -e $infile){
       print "findMaxDistInDiffPDB.pl -p1 $A -p2 $B $extraarg -cutoff $cutoff -outf $infile \n";
       system("findMaxDistInDiffPDB.pl -p1 $A -p2 $B $extraarg  -cutoff $cutoff -outf $infile");
	}
	my $cnt = 0 ;


	my $origa = $A ;
	my $origb = $B ;
	if($unmappedisDNA){
	    my $IsDNAA = IsDNA($A);
	    my $IsDNAB = IsDNA($B);
		die "$A is DNA??? " if(exists $namemap->{$A} && $namemap->{$A} ne "DNA" && $IsDNAA);

		if(exists $namemap->{$B} && $namemap->{$B} ne "DNA" && $IsDNAB ){
		    die  "$B is DNA???  " ;
		}

	    $namemap->{$A} = "DNA" if($IsDNAA && !exists $namemap->{$A});
	    $namemap->{$B} = "DNA" if($IsDNAB && !exists $namemap->{$B});
	}
	$A = $namemap->{$A} if(exists $namemap->{$A});
	$B = $namemap->{$B} if(exists $namemap->{$B});

	return if($A =~ /DNA/i && $B =~ /DNA/i);

	my $forall_commontable = {};

	my @HEL1 ;
	my @HEL2 ;
	my $first = 1 ;

	my $RES1 = {};
	my $RES2 = {};
	my $ifh = util_read($infile);
	my $PDB1 ;
	my $PDB2 ;
	while(<$ifh>){
		## sample line
		##  4R8PI 4R8PL    DT/1/OP2 ARG/325/NH1     2.658       X       X    0
		last if($cnt eq $cntmax);

		#print " JJJ $_ \n";
		my @l = split ;
		my $N = @l ;
		my $distance = $l[4];
		die if(!defined $distance);


	    print $ofh " $A $B ===\n" if($first);
	    $first = 0 ;

		my ($p1,$p2,$a1,$a2,$d,$h1,$h2,$dir) = @l ;
		if(! defined $PDB1){
			$PDB1 = $p1 ;
			$PDB2 = $p2 ;
		}
		else{
			die  "$infile" if($p1 ne $PDB1);
			die "$infile" if($p2 ne $PDB2);

		}
		my ($res1,$num1,$type1) = split "/", $a1 ;
		my ($res2,$num2,$type2) = split "/", $a2 ;

		if(! defined $num2){
			 my $ofherr = util_open_or_append("errfile");
		     print $ofherr "$infile\n";
		}
		$RES1->{$num1} = $res1 ;
		$RES2->{$num2} = $res2;
			
		push @HEL1,$h1  if($h1 ne "X");
		push @HEL2,$h2  if($h2 ne "X");

		if($pairwise){
		    $cnt++;
		    print $ofh $_ ;
		}
		else{
			
			$info->{$p2} = {} if(! defined $info->{$p2});
			$info->{$p2}->{$num1} = {} if(!defined $info->{$p2}->{$num1}) ;
			$info->{$p2}->{$num1}->{$a1} = 1;
			$forall_commontable->{$num1} = "$a1 $a2" ;

		}
		last if($distance > $cutoff);
	}


    my $RESSTR1 = "";
	foreach my $k (sort {$a <=> $b} keys %{$RES1}){
		my $res = $RES1->{$k} ;
		$RESSTR1 = $RESSTR1 . $res . $k . ",";
	}
    my $RESSTR2 = "";
	foreach my $k (sort {$a <=> $b} keys %{$RES2}){
		my $res = $RES2->{$k} ;
		$RESSTR2 = $RESSTR2 . $res . $k . ",";
	}


	@HEL1 = util_sortuniqArray(@HEL1);
	@HEL2 = util_sortuniqArray(@HEL2);

	my $H1str = join ",", @HEL1;
	my $H2str = join",",  @HEL2;
	$H1str = "X" if($H1str eq "");
	$H2str = "X" if($H2str eq "");

	if(!$first){
	    print $ofh "EDGE $PDB1 $PDB2 $A $B $H1str $H2str $RESSTR1 $RESSTR2 ===\n\n\n" ;
		system("mkdir -p PYMOLCOMMANDS");
		if($pymolcommands){
		     my $concat = "$PDB1.$PDB2";
		     my $concatfile = "PYMOLCOMMANDS/$concat.csh";
		     my $OFHPymol = util_write("$concatfile");
		     print $OFHPymol "pymol.2proteins.pl -outf $concat.p1m -pdb1 $PDBDIR/$PDB1.pdb -pdb2 $PDBDIR/$PDB2.pdb\n";
		     print $OFHPymol "pymol.listofresidues.pl -out $concat.p1m -strin $RESSTR1 -chain A -col blue -three -PDB PDBA\n";
		     print $OFHPymol "pymol.listofresidues.pl -out $concat.p1m -strin $RESSTR2 -chain A -col green -three -PDB PDBB\n";
		}
	}

	push @FORFINDINGCOMMON,$forall_commontable;
	close($ifh);
	unlink $infile if($removefile);
}


sub IsDNA{
	my ($protein) = @_ ;
	my $pdb = "$PDBDIR/$protein.pdb";
	my $pdb1 = new PDB();
	$pdb1->ReadPDB($pdb);
	my $IsDNA = $pdb1->IsDNA();
	return $IsDNA ;
}


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
