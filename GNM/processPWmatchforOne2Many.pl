#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use MyGeom;
use MyGNM;
use PDB;
use ConfigPDB;
use Math::Geometry ;
use Math::Geometry::Planar;


use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($genomefile,$infile,$f1,$f2,$outfile,$cutoff,$getcommands,$which_tech,$listfile,$protein);
my ($isNT,$outdir,$ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 0 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "f1=s"=>\$f1 ,
            "f2=s"=>\$f2 ,
            "listfile=s"=>\$listfile ,
            "genomefile=s"=>\$genomefile ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "outdir=s"=>\$outdir ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "isNT=i"=>\$isNT ,
            "getcommands=i"=>\$getcommands ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
my $ofh = util_write($outfile);
my $ofhsinglecommands = util_write("$outfile.single.csh");

usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a listfile -option -f1  ") if(!defined $f1);
usage( "Need to give a listfile -option -f2  ") if(!defined $f2);
usage( "Need to give a listfile -option -isNT  ") if(!defined $isNT);
usage( "Need to give a outdir pdb id -option -outdir  ") if(!defined $outdir);
usage( "Need to give a getcommands ") if(!defined $getcommands);

## genomefile is to get the annotation
usage( "Need to give a output file name => option -genomefile ") if(!defined $genomefile && !$getcommands);
my $map2name = {};
$map2name = util_mapID2fullStringInFastaFile($genomefile) if(!$getcommands);

my $CNT = 0 ; 
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();


my $BLASTCOMMAND = $isNT ? "blastn" : "blastp";
my $ifh = util_read($infile);


my $TARGETS = {};
while(<$ifh>){
	my ($i,@l) = split ;

	my $maxscore = 0 ;
	my $bestmatch = 0 ;
	my $cnt = 0 ;

	my $f2files = "";
	foreach my $j (@l){
		 $TARGETS->{$j} = 1;
		 next if($i eq $j); ## Can happen for self ...
		 $f2files = $f2files . " -expr $f2/$j.ALL.1.fasta ";
		 my $outfile = "$outdir/$i.$j.blast.nt";
		 if($getcommands){
		    if(! -e $outfile){
				$cnt++;
	            print $ofh "$BLASTCOMMAND -query $f1/$i.ALL.1.fasta -subject $f2/$j.ALL.1.fasta -out $outfile\n";
			}
		 }
		 else{

		    my ($info,$querylength,$Subjectname,$queryname,$blastscore,$expect) = GNM_PARSEBLAST_BESTVAL($outfile);
			die "- querylength not defined in $outfile" if(!defined $querylength);
			next if(!defined $blastscore);
			if($blastscore > $maxscore){
			   $maxscore = $blastscore ;
			   $bestmatch = $j ;
			}
			
		}

	}


	my $blastoutfile = "$outdir/$i.blast.nt";
	if(! -e $blastoutfile){
	      print $ofhsinglecommands "processPWmatchforOne2ManySingle.pl -f1 $f1/$i.ALL.1.fasta -tag $i $f2files -isNT $isNT -outfile $blastoutfile \n";
	}

	if(!$getcommands){
		 if($maxscore){
		 	   my $NM = $map2name->{$bestmatch};
			   ## The full name includes the first >
			   $NM =~ s/.//;
	           print $ofh "$i $NM $maxscore \n";
		 }
	}
}

if($getcommands){
    system("scheduleprocessInsertingsleep.pl -inf $outfile -inter 200 -slee 1 -out $outfile.SCH");
}
else{
	system("sortOnLast.pl -in $outfile -rever" );
}

my $junk = util_writelist2file("list.targets",(keys %${TARGETS}));

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
