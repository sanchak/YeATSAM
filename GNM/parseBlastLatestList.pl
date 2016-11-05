#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use MyGeom;
use PDB;
use ConfigPDB;
use MyGNM;
use Math::Geometry ;
use Math::Geometry::Planar;

## disable perl's warning mechanism
no warnings 'recursion';

use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($isNT,$findcharstring,$forWGS,$checkforsame,$infile,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions,$trs,$blastdir);
my $strict = 1 ;
my $verbose = 0 ;

my $blastcutoff ;
my $percentlength ;
my $percentmatched ;
my $percentidentity ;
my $expectlimit ;
my $mapfile ;

my $postfix = "blast.nt";
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "mapfile=s"=>\$mapfile ,
            "protein=s"=>\$protein ,
            "infile=s"=>\$infile ,
            "postfix=s"=>\$postfix ,
            "trs=s"=>\$trs ,
            "checkforsame"=>\$checkforsame ,
            "p1=s"=>\$p1 ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "blastdir=s"=>\$blastdir ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "blastcutoff=i"=>\$blastcutoff ,
            "strict=i"=>\$strict ,
            "forWGS=i"=>\$forWGS ,
            "verbose=i"=>\$verbose ,
            "percentlength=i"=>\$percentlength ,
            "percentmatched=i"=>\$percentmatched ,
            "percentidentity=i"=>\$percentidentity ,
            "findcharstring=i"=>\$findcharstring ,
            "isNT=i"=>\$isNT ,
            "expectlimit=f"=>\$expectlimit ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);


usage( "Need to give a input file name => option -listfile ") if(!defined $listfile);
usage( "Need to give a input file name => option -blastdir ") if(!defined $blastdir);
usage( "Need to give a input file name => option -blastcutoff ") if(!defined $blastcutoff);
usage( "Need to give a input file name => option -forWGS ") if(!defined $forWGS);
usage( "Need to give a input file name => option -findcharstring ") if(!defined $findcharstring);
usage( "Need to give a input file name => option -isNT ") if(!defined $isNT);

if(!$forWGS){
    usage( "Need to give a input file name => option -mapfile ") if(!defined $mapfile);
}

my @list= util_read_list_words($listfile);
my $list = {};
map { s/\s*//g ; $list->{$_} = 1 ; } @list ;

my $N = @list ;

my $HACKDB = "LISTS";
#my $HACKDB = "./";
system ("mkdir -p $HACKDB/$outfile");

my $fname = "$HACKDB/$outfile/$outfile.$blastcutoff" ;
print "Info: All in dir $fname\n";
my $fnameannot = $fname . ".anno";

my $ofh = util_write($fname);
my $ofhanno = util_write($fnameannot);

my $ofhSCAFF2TRS = util_write("scaff2trs.$blastcutoff");
my $ofhCommands = util_write("pairwiseTRS2scaffold.csh");
my $found = 0 ;
my $MAPSscafftoTRS = {};
print STDERR "Info: Cutoff is $blastcutoff=blastcutoff\n";
if($strict){
	my $ofhtmp = util_write("errfiles");
	my $error = 0 ;
    foreach my $i (@list){
        $trs = $i ;
        $infile = "$blastdir/$trs.$postfix";
	    #die "$infile does not exist" if(! -e $infile);
		if( -z $infile || ! -e $infile){
	         print $ofhtmp "$i $infile is empty or nothere\n" ;
			 $error = 1;
		}
    }
	if($error){
		system("wc -l errfiles");
	    die " missing files see errfiles" if($error);
	}
	my $NNN = @list ;
	print "Info: Good, all $NNN files in $listfile exist\n";
}

my $CCC = 0 ;
foreach my $i (@list){
    $trs = $i ;
    $infile = "$blastdir/$trs.$postfix";
	next if(! -e $infile);
	$found++;
    my ($info,$querylength,$Subjectname,$queryname) = GNM_PARSEBLAST($infile,$verbose);
	next if($strict eq 0 && !defined $querylength);
    die "$trs - querylength not defined in $infile" if(!defined $querylength);
    my $x = GNM_MakeGroupOfNR_OrAnnotate($info,$querylength,$Subjectname,$isNT,$findcharstring,$forWGS,$trs,$infile,$ofh,$ofhanno,$blastcutoff,$MAPSscafftoTRS,$checkforsame);
}

print "The forWGS arg controls grouping...\n";

print "There are $N items of which found $found \n";
if(!$found){
	print "Warning: Check the postfix .... now it is $postfix\n";
}

if(defined $mapfile){
    system ("annFinalCommands.pl -listfile $listfile -forWGS $forWGS -outf $outfile -blastcutoff $blastcutoff -mapfile $mapfile");
}
else{
    system ("annFinalCommands.pl -listfile $listfile -forWGS $forWGS -outf $outfile -blastcutoff $blastcutoff");
}


$, = "  ";
foreach my $scaff (keys %{$MAPSscafftoTRS}){
	my @trs = (keys %{$MAPSscafftoTRS->{$scaff}}); 
	print $ofhSCAFF2TRS "parseBlastLatestpairwiseAllScaffolds.pl $scaff @trs \n";
}



sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
