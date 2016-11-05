#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use MyPymol;

use PDB;
use Atom;
use Residue;

use POSIX qw(floor);
use Math::Combinatorics;
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($pdb1,$pdb2,$infile,$listfile,$outfile,$dontrunPymol,$dontshowsurface);
my ($name1,$name2);
my @pdbs ;
my @distances ;
GetOptions(
            "pdb1=s"=>\$pdb1 ,
            "pdb2=s"=>\$pdb2 ,
            "listfile=s"=>\$listfile ,
            "dontrunPymol"=>\$dontrunPymol ,
            "dontshowsurface"=>\$dontshowsurface ,
            "infile=s"=>\$infile ,
            "pdb=s"=>\@pdbs,
            "dist:s"=>\@distances,
            "outfile=s"=>\$outfile,
            "name1=s"=>\$name1 
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a output file name => option -listfile ") if(!defined $listfile);
die "Please give .p1m postfix to outfile" if(!($outfile =~ /\.p1m/));
my $ofh = util_write($outfile);
print "Info: Writing $outfile\n";

my @list= util_read_list_sentences($listfile);
my ($RESULTDIR,$PDBDIR) = util_SetEnvVars();

my $ofhcommands = util_write("commands.charge");


my $colorcoding = {};
$colorcoding->{"R"} = "blue";
$colorcoding->{"K"} = "blue";
$colorcoding->{"H"} = "blue";
$colorcoding->{"D"} = "red";
$colorcoding->{"E"} = "red";
my $cnt = 0 ; 
foreach my $i (@list){
	my $pdbfle = "$PDBDIR/$i.pdb" ;
	my $pdb = new PDB($pdbfle);
	$cnt++;
	my $ID = "PROT$cnt";

	my $tmpfh = util_write("tmp");
			 
	my $seq = "";
	my $CNTres = 0 ;
	foreach my $residue (@{$pdb->{RESIDUES}}){
		my $x  = $residue->PrintSingleLetter($pdb);
		next if($x eq "?");
		if(defined $x && exists $colorcoding->{$x}){
			 #$seq = $seq . $x ;
			 my $idx = $residue->GetIdx();
			 my $color = $colorcoding->{$x} ;
			 next if ($color ne "blue");
			 $CNTres++;
			 print $ofhcommands "select block_query$CNTres, /$ID//A/$idx\n";
			 print $ofhcommands "color $color, block_query$CNTres\n";
			 print $ofhcommands "show spheres, block_query$CNTres\n";
		}
	}

	push @pdbs,$pdbfle ;
}

util_PrintPymolManyProteins($ofh,@pdbs);
system("wc -l commands.charge");


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
