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


use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($organism,$idx,$infile,$listfile,$postfix,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@listfiles);
my $howmany = 100000 ;
my $verbose = 0 ;
GetOptions(
            "organism=s"=>\$organism ,
            "protein=s"=>\$protein ,
            "postfix=s"=>\$postfix ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "listfiles=s"=>\@listfiles,
            "howmany=i"=>\$howmany ,
            "idx=i"=>\$idx ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
my $PWD = cwd;
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a input file name => option -organism ") if(!defined $organism);
#my $ifh = util_read($infile);
#usage( "Need to give a input file name => option -fastafile ") if(!defined $fastafile);
#my ($tablefasta,$infoSeq2PDB,$mapChainedName2Name) = util_parsePDBSEQRESNEW($fastafile,0,0);

 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();

my $MINCOUNT = 0 ;
my $MINLENGHT = 100 ;
my $MAXAADIFF = 45 ;
my $CNT = 0 ;

my $countStr = Config_getCountStr();
# for reading fasta
#my ($fe,$x) =  util_readfasta($infile);
#chomp $x ;


my ($tableMETA) = util_maketablefromfile("list.meta");
$CNT = (keys %{$tableMETA});
print "There are $CNT transcripts \n";


#save the complete count as a table with the TRS
my ($tableCount) = util_maketablefromfile_firstentry("COUNTS/bwa_counts_run1.txt.0");
print "There are $CNT trs with cnt > $MINCOUNT \n";

my ($tableORGANISMS) = util_maketablefromfile($organism);

my $ofh = util_write($outfile);

my $sort = {};
foreach my $OO (sort keys %{$tableORGANISMS}){
    my $tableTRSanno = GetTRSwithString("MAPDIRS/map.checknovel.150.anno.real",$OO);
	my (@l) = (keys %{$tableTRSanno});
	my $NUMOFTRS = @l ;
    my $countStr = Config_getCountStr();
    
    my $COUNTS = {};
    foreach my $i (@l){
    	my $s = $tableCount->{$i} ;
        my (@lvalues) = split " ",$s;
        my (@lstrs) = split " ",$countStr;
		my $NM = $lvalues[0];
		my $SSS = $tableTRSanno->{$NM} ;
		print $ofh "$OO $s\n" if($verbose);
		print $ofh "$OO $SSS\n" if($verbose);
    
    	shift @lvalues ;
    	shift @lstrs ;
    	while(@lstrs){
    		my $id = shift @lstrs ;
    		my $v = shift @lvalues ;
    		if(! exists $COUNTS->{$id}){
    			$COUNTS->{$id} = 0 ;
    		}
    		$COUNTS->{$id} = $COUNTS->{$id} + $v ;
    	}
    }
    
    
	my $finalcnt = 0 ; 
	my $finalstr = "$OO NUMOFTRS=$NUMOFTRS-----------\n";
    foreach my $id (sort keys %{$COUNTS}){
       my $v = $COUNTS->{$id} ;
       #if($v){
       	 $finalstr = $finalstr .  "($id, $v) \n";
       #}
	   $finalcnt = $finalcnt + $v ;
    }
	$sort->{$finalstr} = $finalcnt ;
}

my @sorted = sort { $sort->{$b} <=> $sort->{$a}} (keys %{$sort});
foreach my $i (@sorted){
	print $ofh $i ;
}

sub GetTRSwithString{
	my ($in,$str) = @_ ; 
	my $ifh = util_read($in);
	my $table = {};
	while(<$ifh>){
		if(/$str/i){
			my ($nm) =split ;
			$table->{$nm} = $_ ;
		}
	}
	return $table ;
}


###############
sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
