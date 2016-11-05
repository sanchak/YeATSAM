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
my ($ispolar,$infile,$outfile,$initial,$which_tech,$hetatm,$listfile,$protein);
my (@expressions);
my $howmany = 4 ;
my $maxdist = 10 ;
my $verbose = 0 ;
my $INCR = 0.1 ;

my $IGNOREBACKBONE = 0 ;
my $doall = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "doall=i"=>\$doall ,
            "listfile=s"=>\$listfile ,
            "hetatm=s"=>\$hetatm ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "maxdist=i"=>\$maxdist ,
            "ispolar=i"=>\$ispolar ,
            "initial=f"=>\$initial ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
## for debug 
usage( "Need to give a hetatm -option -hetatm  ") if(!defined $hetatm);
usage( "Need to give a protein pdb id -option -protein  ") if(!defined $protein);
usage( "Need to give a ispolar pdb id -option -ispolar  ") if(!defined $ispolar);
my $fhallinfo = util_write("$protein.allinfo");
my $CNT = 0 ; 
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();
my $PWD = cwd;

my $DONEATOMS = {};

my $pdbfile = "$PDBDIR/$protein.pdb";
my $pdb1 = new PDB($pdbfile);
$pdb1->ReadPDB($pdbfile);



my $infohetsTOP = $pdb1->FindHetAtms();
die "Error: No hetatms for $hetatm found in $pdbfile" if(!exists $infohetsTOP->{$hetatm});
my $infohets = $infohetsTOP->{$hetatm};

my $Nhets = (keys %{$infohets});
my $OUTDIR = "HETINFO/$protein/$hetatm";
system("mkdir -p $OUTDIR");


## Hydrogen bonds
my $Allowed = {};
$Allowed->{"ON"} = 1 ; $Allowed->{"NO"} = 1 ;
$Allowed->{"OH"} = 1 ; $Allowed->{"HO"} = 1 ;
$Allowed->{"OO"} = 1 ;
$Allowed->{"NN"} = 1 ;
$Allowed->{"NH"} = 1 ; $Allowed->{"HN"} = 1 ;
$Allowed->{"SH"} = 1 ; $Allowed->{"HS"} = 1 ;

my @backboneatoms = qw( O N CA);
my $backboneatomstable = util_make_table(\@backboneatoms);

my @NONPOLAR = qw(PRO ALA GLY ILE LEU MET VAL);
my $NONPOLAR = {};
foreach my $i (@NONPOLAR){
	$NONPOLAR->{$i}= 1;
}



my $log = "$OUTDIR/$protein.log";
my $ofhlog = util_write($log);

my $pseudoatom = new Atom();
$pseudoatom->SetIdx(10000);


my $RTouched = {} ;
my $RTouchedResAtom = {} ;
my $RTouchedDrugAtom = {} ;


## There will be multiple lines here
my $claspout = "$OUTDIR/$protein.$howmany.clasp.in";
my $table = "$OUTDIR/$protein.$howmany.$ispolar.table.in";
my $fhtable = util_write($table);
my $fhclasp = util_write($claspout);


print "There are $Nhets instance of $hetatm\n";
foreach my $resnum (sort { $a <=> $b} keys %{$infohets}){
	my @atomlines = @{$infohets->{$resnum}};
	my $Natoms  = @atomlines ;
	print "Info: Resnum $resnum for hetatm $hetatm has $Natoms \n";
    my $newsvrpdb = "$OUTDIR/$protein.$hetatm.$resnum.pdb";
    my $ofhpdb = util_write($newsvrpdb);
	foreach my $atomline (@atomlines){
        my($atomstr,$serialnum,$atomnm,$alt_loc,$resname,$chainId,$resnum,$codeforinsertion,@coords)=util_ReadLine($atomline);
	    print $ofhpdb "$atomline \n"; ;
	}

    my $R = defined $initial ? $initial : 2.5 ;
    my ($N);
    my (@t);
    while((keys %{$RTouched}) < $howmany){
	    $R = $R + $INCR ;
	    $N = keys %{$RTouched} ;
	    @t = (keys %{$RTouched});
	    print "rad = $R number found: $N Res numbers: @t ============ \n";
	    ProcessOneRadiiforHetatm($R,$RTouched,$RTouchedResAtom,$RTouchedDrugAtom,$hetatm,@atomlines);
	    last if($R > $maxdist);
    }



    $N = keys %{$RTouched} ;
    @t = (keys %{$RTouched});
    print "rad = $R number found: $N Res numbers: @t ============ \n";
    print "Found $N residues in $R radius\n";



    my $cnt = 0 ;
    foreach my $key (sort {$RTouched->{$a} <=> $RTouched->{$b}} keys %{$RTouched}){
	       $cnt++; 
           last if($cnt > $howmany);
	       my ($res) = $pdb1->GetResidueIdx($key);
	       my $v = $RTouched->{$key};
	       $v =~ s/\s//g;
	       my $atomnm = $RTouchedDrugAtom->{$key};
	       my $type = $RTouchedResAtom->{$key};
	       my $nm = $res->GetName() . "$key";;
	       print $fhclasp "$nm ";
	       print $fhtable "$nm/$type/$atomnm/$v ";
    }
    print $fhclasp "\n";
    print $fhtable "\n";
    print "Wrote in $claspout \n";
    
}

##################################################




sub ProcessOneRadiiforHetatm{
   my ($radii,$residuestouched,$residuestouchedResatom,$residuestouchedDrugatom,$hetatmnm,@atomlines) = @_ ;
   
   my $info = {};
   foreach my $atomline (@atomlines){
     my($atomstr,$serialnum,$atomnm,$alt_loc,$resname,$chainId,$resnum,$codeforinsertion,$x,$y,$z)=util_ReadLine($atomline);

	 $pseudoatom->SetCoords($x,$y,$z);

     my $list = util_make_list($pseudoatom);
	 my ($junk,$neigh)  = $pdb1->GetNeighbourHoodAtom($list,$radii);



	 my $sort = {} ;
	 my $sortedneigh = {};
	 foreach my $a (@{$neigh}){
		my $d = $pdb1->DistanceAtoms($pseudoatom,$a);
	 	my $nm = $a->GetName();
		$sortedneigh->{$a} = $d ;
	 }
	 my @sorted = sort { $sortedneigh->{$a} <=> $sortedneigh->{$b} } (@{$neigh});
	 foreach my $a (@sorted){

			   my $NNMM = $a->GetName();
		       my $d = $pdb1->DistanceAtoms($pseudoatom,$a);

			   print "$NNMM $atomnm $d lllllllllll\n" if($verbose);

			   my $num = $a->GetResNum();
			   my $type = $a->GetType();
	           my ($res) = $pdb1->GetResidueIdx($num);
			   my $resnm = $res->GetName();


	           next if($IGNOREBACKBONE && exists $backboneatomstable->{$type});


			   my @typel = split "", $type ;
			   $atomnm =~ s/\s//g;
			   my ($atomfirst) = ($atomnm =~ /(.)/);
			   my ($typefirst) = ($type =~ /(.)/);
			   my $X = $atomfirst . $typefirst ;
				
			   ## for debug 
		   	   if(! exists $DONEATOMS->{$a}){
			       print $fhallinfo "$resnm$num$type $atomnm $X $d\n";
			   }
			   $DONEATOMS->{$a} = 1 ;

			   if(!$doall){
			       if($ispolar){
			           next if(exists $NONPOLAR->{$resnm});
			           next if(!exists $Allowed->{$X});
			       }
			       else{
			           #next if(!exists $NONPOLAR->{$resnm});
			           next if(exists $Allowed->{$X});
			       }
			   }

	                    
			   if(! exists $residuestouched->{$num}){
			      $residuestouched->{$num} = util_format_float($d,1) ;
			      $residuestouchedResatom->{$num} = $type ;
			      $residuestouchedDrugatom->{$num} = $atomnm ;
			   }
			   my $nm = $a->GetName();
			   $a->Print() if($verbose);
               my $datafile = "$OUTDIR/$hetatmnm.$radii.dat";
               my $ofhdata ;
               if($verbose){
                    $ofhdata = util_write($datafile);
			        print $ofhdata "$atomnm $protein $resnm $num $type $d \n";
					close($ofhdata);
               }
			   print $ofhlog "$radii $atomnm $protein $resnm $num $type $d \n";
     } ## sorted printing
  }  # foreach hetam line

}




sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}

