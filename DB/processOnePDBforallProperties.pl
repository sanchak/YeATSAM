#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
  use Time::HiRes qw( usleep ualarm gettimeofday tv_interval
   clock_gettime clock_getres  clock
   );


use PDB;
use Atom;
use Residue;

use POSIX qw(floor);
use Math::Combinatorics;
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($pdb1,$pdb2,$outfile,$atomidx,$dontrunpymol);
my ($interactive,$annotate,$dist,$mutatefile,$maxresults,$inconf,$outconf,$resultfile,$checkself);
my ($grpconfig) = $ENV{CONFIGGRP} or die ;
my ($justseqeunce,$aalist,$cutoff,$exception,$ann,$config,$p1,$p2,$infile,$score,$ignorepro,$which_tech,$listfile,$protein);
my $MINDIST = 2 ;
my $force = 0 ;
my $clean = 0 ;
my $writeIndivual = 0 ;
$, = "  ";
GetOptions(
            "protein=s"=>\$protein ,
            "mutatefile=s"=>\$mutatefile ,
            "which_tech=s"=>\$which_tech ,
            "p1=s"=>\$p1 ,
            "aalist=s"=>\$aalist ,
            "p2=s"=>\$p2 ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "ann=s"=>\$ann ,
            "cutoff=f"=>\$cutoff ,
            "config=s"=>\$config,
            "force=i"=>\$force,
            "clean=i"=>\$clean,
            "writeIndivual=i"=>\$writeIndivual,
            "justseqeunce=s"=>\$justseqeunce,
            "score=s"=>\$score,
            "exception=s"=>\$exception ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);

my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP) = util_SetEnvVars();


usage( "Need to give protein => option -protein ") if(!defined $protein);
usage( "Need to give protein => option -outfile ") if(!defined $outfile);
if($clean){
	unlink $outfile ;
}

my $origpdb = $protein ;
my $pdb = "$PDBDIR/$protein.pdb";

$pdb1 = new PDB();
$pdb1->ReadPDB($pdb);


my $fastafile = "$FASTADIR/$origpdb.ALL.1.fasta";
my $ofh = util_write($fastafile);
my $seq = $pdb1->WriteFasta($protein,$ofh);
my $seqlen = length($seq);
if(!$seqlen){
	die "$protein has no residues\n";
}
print STDERR "Wrote fasta in $fastafile of length $seqlen\n";


## Now helices - get PDBs from DSSP
my $dsspinfile = "$DSSP/$protein.dssp";
die "Please run DSSP " if(! -e $dsspinfile);

# We do not need beta processing for the time being...
#util_ParseDSSP($protein,$pdb1,$dsspinfile,"BETA","finaloutfile",$writeIndivual,$force);

my $helixlist = util_ParseDSSP($protein,$pdb1,$dsspinfile,"HELIX",$outfile,$writeIndivual,$force);



my $cntdisulfhide = util_GetDisulphide ($protein,$pdb1,1);
print "There were $cntdisulfhide disulphide bonds\n";

$pdb1->ProcessDisulphideForHelices();


print STDERR "Info: writeIndivual = $writeIndivual\n";
if($writeIndivual){
    SortHelixValues($protein,"POS");
    SortHelixValues($protein,"NEG");
    SortHelixValues($protein,"HYD");
    $pdb1->PrintABSequence();
}

if(0){
system("cd $PDBDIR; fpocket -f $protein.pdb  ; cd -");
my ($infofpocket,$infoperatom,$volume) = parseFPocketFile($protein,$pdb1,1,"$PDBDIR/${protein}_out/pockets/pocket0_atm.pdb");
util_table_print($infofpocket);
($infofpocket,$infoperatom,$volume) = parseFPocketFile($protein,$pdb1,1,"$PDBDIR/${protein}_out/pockets/pocket1_atm.pdb");
util_table_print($infofpocket);
system("rm -rf $PDBDIR/${protein}_out/");
}




sub SortHelixValues{
   my ($protein,$whattodo) =@_;
   system("helixSortResults.pl -outf $protein.$whattodo -inf $HELIXDIR//$protein.HELIXVALUES -what $whattodo");
   system("echo -n \"$whattodo \" ; head -1 $protein.$whattodo");
   system("echo ==========");
}


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
