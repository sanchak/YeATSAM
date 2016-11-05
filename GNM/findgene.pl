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
use AAConfig;



use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($increment,$dbfile,$idx,$infile,$postfix,$p1,$p2,$outfile,$cutoff,$tag,$listfile,$blastdir,$protein);
my ($ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "tag=s"=>\$tag ,
            "blastdir=s"=>\$blastdir ,
            "protein=s"=>\$protein ,
            "dbfile=s"=>\$dbfile ,
            "postfix=s"=>\$postfix ,
            "infile=s"=>\$infile ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "idx=i"=>\$idx ,
            "increment=i"=>\$increment ,
            "verbose=i"=>\$verbose ,
            "cutoff=i"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a input file name => option -dbfile ") if(!defined $dbfile);
usage( "Need to give a input file name => option -tag ") if(!defined $tag);

usage( "Need to give a protein pdb id -option -protein  ") if(!defined $protein);
usage( "Need to give a cutoff pdb id -option -cutoff  ") if(!defined $cutoff);
usage( "Need to give a blastdir pdb id -option -blastdir  ") if(!defined $blastdir);
usage( "Need to give a increment pdb id -option -increment  ") if(!defined $increment);
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();

system ("mkdir -p FINDGENETMP");

my $ofh = util_write($outfile);
system ("echo $protein >  $outfile");
close($ofh);


my $prev ;
my $CNT = 0 ;
system("cp -f $outfile FINDGENETMP/addedInStep.$CNT");
while(1){
	$CNT++;
    my @list= util_read_list_sentences($outfile);
    my $N = @list;
    print "There are $N blasts to run\n";
    

    my $ofhcomm = util_write("$tag.commands.csh");
	print $ofhcomm "#!/bin/csh -f\n";

    foreach my $i (@list){
       if(! -e "$blastdir/$i.blast.nt"){
          print $ofhcomm "BLASTP $dbfile $FASTADIR/$i.ALL.1.fasta $blastdir/$i.blast.nt \n";
       }
    }
    close($ofhcomm);
    system ("wc -l $tag.commands.csh");
    system ("cat $tag.commands.csh");
    system ("chmod +x $tag.commands.csh");
    system ("./$tag.commands.csh");

    system("parseBlastLatestList.pl -out FINDGENETMP/$tag -lis $outfile -blastdir $blastdir/ -blastcutoff $cutoff -forWGS 0 -findcha 0 -isNT 0 -str 1");
    system("printInColumnUniquely.pl -inf FINDGENETMP/$tag.$cutoff.appended.removevalues");
    system("twolists.pl FINDGENETMP/$tag.$cutoff.appended.removevalues.singleline $outfile");
    system("cp -f ofhinAbutnotinB FINDGENETMP/addedInStep.$CNT");
    system("cp -f FINDGENETMP/$tag.$cutoff.appended.removevalues.singleline $outfile");
	$cutoff = $cutoff + $increment ;
    my @list= util_read_list_sentences($outfile);
	my $NEWN = @list ;
	if($NEWN < $N ){
		print "Info:Finishing off at $NEWN since equal to prev $N, at cutoff $cutoff\n";
		last ;
	}

	#die if ($CNT eq 2);
}



sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
