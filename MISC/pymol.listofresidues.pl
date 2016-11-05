#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use ConfigPDB;
use MyGeom;
use PDB;
use MyConfigs;

use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($stringresidue,$threeCodeRemove,$singleCodeRemove,$PDBID,$chainID,$infile,$outfile,$onecolor,$color,$listfile,$protein);
my (@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "color=s"=>\$color ,
            "chainID=s"=>\$chainID ,
            "PDBID=s"=>\$PDBID ,
            "stringresidue=s"=>\$stringresidue ,
            "infile=s"=>\$infile ,
            "listfile=s"=>\$listfile ,
            "outfile=s"=>\$outfile ,
            "singleCodeRemove"=>\$singleCodeRemove ,
            "threeCodeRemove"=>\$threeCodeRemove ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
die "Expect a .p1m file $outfile to append to"  if(! -e $outfile);
my $ofh = util_append($outfile);
usage( "Need to give a listfile or stringresidue -option -listfile / -stringresidue  ") if(!defined $listfile && !defined $stringresidue);
usage( "Need to give one of listfile or stringresidue -option -listfile / -stringresidue  ") if(defined $listfile && defined $stringresidue);
usage( "Need to give a PDBA or PDBB pdb id -option -PDBID  ") if(!defined $PDBID);
usage( "Need to give a chainID pdb id -option -chainID  ") if(!defined $chainID);
usage( "Need to give a color pdb id -option -color  ") if(!defined $color);
my $what = defined $listfile? $listfile: $stringresidue;
print $ofh "\n\n# Info: Appending to $outfile from $what for chainID $chainID, PDBID $PDBID\n";
print $ofh "\n\nset sphere_scale, 0.60, (all)\n";

my $CNT = 0 ; 
my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
my $PWD = cwd;

my @colors = qw ( magenta cyan yellow blue black red green gray orange  ) ;


print "Info: Color coding is - Red=+ve, magenta=-ve. Shape - charged is spheres, else sticks.\n";

## lsit contains list of 2 numbers - start and end
my @list ;
my $junktable ;
if(defined $listfile){
   @list = util_read_list_words($listfile);
}
else{
   ($junktable,@list) = util_CreateTableAndListFromSymbolSeperated($stringresidue,",");
}
my $list = {};
map { s/\s*//g ; $list->{$_} = 1 ; } @list ;
@list = (keys %{$list});
my $N = @list ; 
print "Info: Appending $N residues to $outfile\n";

my $cnt = 0 ; 
my $origcolor = $color ;
while(@list){

    ## start with the orig color, it changes in the code - so need to restore
	$color = $origcolor ;

   my $a = shift @list ;

   my ($tableATOM,$HYDROVAL,$colortable,$value,$chargedtable,$aromatic) = Config_Helix();
   my $charged = 0 ;
   if(defined $singleCodeRemove){
	   my ($res) = ($a =~ /(.)/);
	   if(exists $chargedtable->{$res}){
	   	  $charged = 1;
		  $color = $chargedtable->{$res} eq 1 ? "red" : "magenta";
	   }
	   elsif (exists $aromatic->{$res}){
	   	  $color = "yellow";
	   }
       $a =~ s/.// ;
   }

   if(defined $threeCodeRemove){
	   my ($res) = ($a =~ /(...)/);
	   if(exists $chargedtable->{$res}){
	      $charged = 1 ;
		  $color = $chargedtable->{$res} eq 1 ? "red" : "magenta";
	   }
	   elsif (exists $aromatic->{$res}){
	   	  $color = "yellow";
	   }
       $a =~ s/^...// ;
   }

   my $whattodraw =  $charged ? "spheres" : "sticks";
   print $ofh "select block_query$cnt, /$PDBID//$chainID/$a\n";
   print $ofh "color $color, block_query$cnt\n";
   print $ofh "show $whattodraw, block_query$cnt \n";

}


chmod 0777, $outfile ;

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
