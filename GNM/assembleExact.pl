#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($all,$infile,$outfile,$or,$silent,$groupinfo);
my ($size,$DIR,$listfile,$ignorefile,$mergedir);
my $howmany = 600000 ; 
my $cutofflength = 0 ; 
my @types = (); 
my $protein ;
my @motifs = (); 
GetOptions(
            "all"=>\$all ,
            "groupinfo"=>\$groupinfo ,
            "silent"=>\$silent ,
            "infile=s"=>\$infile ,
            "mergedir=s"=>\$mergedir ,
            "dir=s"=>\$DIR ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "howmany=i"=>\$howmany ,
            "size=i"=>\$size ,
            "or=i"=>\$or ,
            "cutofflength=i"=>\$cutofflength ,
            "type=s"=>\@types,
            "protein=s"=>\$protein,
            "motif=s"=>\@motifs,
            "outfile=s"=>\$outfile 
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
print STDERR "Info: parsing file $infile - might take some time\n";
my ($info,$infoSeq2PDB,$mapChainedName2Name) = util_parsePDBSEQRESNEW($infile,0);

## map seq and names to numbers
my $TRSID = 0 ;
my $MAPCNT2NAME = {};
my $MAPCNT2SEQ = {};
my $MAPSEQ2CNT = {};

print "Mapping names to numbers\n";
foreach my $seq (keys %{$infoSeq2PDB}){

	my @l = sort @{$infoSeq2PDB->{$seq}};
	my $name = $l[0];
	$MAPCNT2SEQ->{$TRSID} = $seq ;
	$MAPCNT2NAME->{$TRSID} = $name ;
	$MAPSEQ2CNT->{$seq} = $TRSID ;
	$TRSID++;
}

my $done = {};
my $total = 0 ;
my $final = 0 ;

my $BEGIN = {};
my $END = {};

### Concat all end strings
my $ENDSTR = "";

my $TABforENDSTR = {};
my $TABforBEGINSTR = {};

my $multipleBEGIN = {};
my $multipleEND = {};
foreach my $seq (keys %{$MAPSEQ2CNT}){
	my $len = length($seq);
	my $TRSID = $MAPSEQ2CNT->{$seq} ;

	## first n and last n 
	my ($begin,$end) = util_GetTerminalStrings($seq,100);
	$BEGIN->{$TRSID} = $begin;
	$END->{$TRSID} = $end;

	if(! defined $TABforENDSTR->{$end} ){
	   $TABforENDSTR->{$end} = [] ;
	}
	else{
		$multipleEND->{$end} = 1 ;
	}
	push @{$TABforENDSTR->{$end}} , $TRSID ;

	if(! defined $TABforBEGINSTR->{$begin} ){
	     $TABforBEGINSTR->{$begin} = [] ;
	}
	else{
		$multipleBEGIN->{$begin} = 1 ;
	}
	push @{$TABforBEGINSTR->{$begin}} , $TRSID ;
}

my $Nbegin = (keys %{$multipleBEGIN});
my $NEnd = (keys %{$multipleEND});

print  "There are Nbegin = $Nbegin and $NEnd =NEnd\n";

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
