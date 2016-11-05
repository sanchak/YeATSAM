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
my ($writefasta,$trs,$orfdir,$pdbseqres,$infile,$outfile,$which_tech,$listfile,$protein);
my (@expressions);
my $size ;
my $verbose = 1 ;
GetOptions(
            "writefasta=i"=>\$writefasta ,
            "pdbseqres=s"=>\$pdbseqres ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "size=i"=>\$size ,
            "verbose=i"=>\$verbose ,
            "orfdir=s"=>\$orfdir ,
            "trs=s"=>\$trs ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);

usage( "Need to give a input file name => option -trs ") if(!defined $trs);
usage( "Need to give a input file name => option -orfdir ") if(!defined $orfdir);


my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT) = util_SetEnvVars();
my $PWD = cwd;


my $ofhlongest = util_open_or_append("list.orflongest");
$pdbseqres = "$orfdir/$trs.orf";
print "READING $pdbseqres\n" if($verbose);
my ($info,$infoSeq2PDB,$mapChainedName2Name) = util_parsePDBSEQRESNEW($pdbseqres,0);

my $ofherror = util_open_or_append("list.err.repeatorf");

my $sort = {};
foreach my $k (keys %{$infoSeq2PDB}){
	my $NEWK = $k ;
	## Remove X while measuring length
	$NEWK =~ s/X//g;
	my $l = length($NEWK);
	$sort->{$k} = $l ;
}

my @sorted = sort { $sort->{$b} <=> $sort->{$a}} (keys %{$sort});
	

system ("mkdir -p FASTADIR_ORF");
system ("mkdir -p FASTADIR_JUSTONLONGEST");
my $done = {};
my $cnt = 0 ;
my $errorprinted = 0 ;
my $ofhlist = util_write("FASTADIR_ORF/$trs.list");
foreach my $k (@sorted){
	die if(exists $done->{$k});
	$done->{$k} = 1 ;
	#print "$k\n";

	$cnt++;
	my $l = length($k);

	my @l = @{$infoSeq2PDB->{$k}};
    my $N = @l ;
	if($N > 1 && !$errorprinted){
		my $LEN = length($k);
		print $ofherror "$trs $cnt n=$N $k $LEN\n";
		$errorprinted = 1;
	}

	### writefasta has the number of ORFS to write 
	if(defined $writefasta && $cnt < $writefasta){
		my @lll = split " ", $l[0];
		my $name = "FASTADIR_ORF/$lll[0]" . ".ALL.1.fasta";
		print $ofhlist "$lll[0]\n";
		my $len = length($k);
		print "Writing $name $len\n" if($verbose);
		if(! -e $name){
		   my $ofh = util_write($name);
	           print $ofh ">$l[0] \n";
	           print $ofh  "$k\n";
		   close($ofh);
		}
		if($cnt eq 1){    
	        print $ofhlongest "$lll[0]\n";
		    my $name = "FASTADIR_JUSTONLONGEST/$trs" . ".ALL.1.fasta";
		    print "\tWriting longest $name  $lll[0] $len\n" if($verbose);
		    if(! -e $name){
		       my $ofh = util_write($name);
	           print $ofh ">$trs $lll[0] \n";
	           print $ofh  "$k\n";
		       close($ofh);
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
