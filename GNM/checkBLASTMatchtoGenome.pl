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
my ($idx,$infile,$postfix,$p1,$p2,$outdir,$cutoff,$fastafile,$which_tech,$listfile,$DB,$outfile);
my ($ignorefile,@expressions);
my $howmany = 100000 ;
my $cutoff = 250 ;
my $savedata = 0 ;
my $verbose = 1 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "DB=s"=>\$DB ,
            "postfix=s"=>\$postfix ,
            "infile=s"=>\$infile ,
            "fastafile=s"=>\$fastafile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "outdir=s"=>\$outdir ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "savedata=i"=>\$savedata ,
            "idx=i"=>\$idx ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outdir ") if(!defined $outdir);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_open_or_append($outfile);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
my $ifh = util_read($infile);
usage( "Need to give a input file name => option -fastafile ") if(!defined $fastafile);
my ($tablefasta,$infoSeq2PDB,$mapChainedName2Name) = util_parsePDBSEQRESNEW($fastafile,0,0);

usage( "Need to give a DB pdb id -option -DB  ") if(!defined $DB);
 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();

my ($ignoretable)  = util_maketablefromfile($outfile);


my @list ; 
if(defined $listfile){
   @list= util_read_list_sentences($listfile);
}
else{
	@list = (keys %{$tablefasta});
}

foreach my $i (@list){
	next if(exists $ignoretable->{$i});
	my $OUTF = "$outdir/$i.blast";

	my $infasta = "$FASTADIR/$i.ALL.1.fasta";
	my $seq = $tablefasta->{$i} ; 
	if(! -e $infasta){
		my $ofhfasta = util_write($infasta);
		print $ofhfasta ">$i\n";
		print $ofhfasta "$seq\n";
	}
	system ("blastn -db $DB -query $infasta -out $OUTF");
	my ($found,$score) = ParseBlastQuick($OUTF,$cutoff);
    print $ofh "$i $score $found\n";

	if(! $savedata){
		unlink ($OUTF);
		unlink ($infasta);
	}
}


sub ParseBlastQuick{
	my ($blastfile,$cutoff) = @_ ;
	my $ifh = util_read($blastfile);

	my $found = 0 ;
	while(<$ifh>){
		if(/Sequences producing significant alignments/){
			$_ = <$ifh> ;
			$_ = <$ifh> ;
			#print "$blastfile $_ \n"; 
			my @l = split ;
			my $N = @l - 2;
			my $score = $l[$N];
		    close($ifh);
			if($score > $cutoff){
			     return (1,$score) ;
			}
			return (0,$score) ;
		}
		elsif(/No hits found/){
		    close($ifh);
			return (0,-1) ;
		}
	}
	die "Should not be here";
	
	
}

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
