#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use MyGeom;
use MyGNM;
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
my ($fastafile,$idx,$infile,$postfix,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
$cutoff = 60 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "protein=s"=>\$protein ,
            "fastafile=s"=>\$fastafile ,
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
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a input file name => option -outfile ") if(!defined $outfile);
my $ofhlist = util_write("$outfile");
my $ofhlistsmall = util_write("$outfile.smaller.$cutoff");
my $ofhlistlongest = util_write("$outfile.longest");

my ($NAME,$sequence);
my $info = {};
my $seenNAMES = {};
my $docheck = 1 ;
## have a first parse that save the length...
my $ifh = util_read($infile);
while(<$ifh>){
	     next if(/^\s*$/);
         if(/^\s*>/){
	 	    if(defined $NAME){
			   if(! exists $seenNAMES->{$NAME} ){
	               $docheck = ProcessSingleSequence($NAME,$sequence,$docheck);
			       $seenNAMES->{$NAME} = 1 ;
			   }
			   else{
			        warn "Warning: repeat name $NAME\n" ;
			   }
               $sequence = "";
		    }
    
		    s/>//;
		    ($NAME) = split ;
	     }
	     else{
	        die if(/>/);
		    ## remove spaces 
		    s/ //g;
		    chomp;
		    $sequence = $sequence . $_ ;
	     }
}
if(! exists $seenNAMES->{$NAME} ){
    $docheck = ProcessSingleSequence($NAME,$sequence,$docheck);
}
close($ifh);



## Process each ID - and choose the top three 
my $best = {};
foreach my $key (sort keys %{$info}){
	#print "$key ...\n";
	my $tab = $info->{$key} ;

	my @sl = sort { $tab->{$b} <=> $tab->{$a} } keys %{$tab} ;
	#print "@sl \n";
	$best->{$key} = {};
	my $bestlength = $tab->{$sl[0]} ;
	if($bestlength < $cutoff){
		print $ofhlistsmall "$key $bestlength\n";
	}
	$best->{$key}->{$sl[0]} = 1;
	$best->{$key}->{$sl[1]} = 2;
	$best->{$key}->{$sl[2]} = 3;
}

$ifh = util_read($infile);
undef $NAME ;
while(<$ifh>){
	     next if(/^\s*$/);
         if(/^\s*>/){
	 	    if(defined $NAME){
	           $docheck = PrintSequence($NAME,$sequence,$docheck);
               $sequence = "";
		    }
    
		    s/>//;
		    ($NAME) = split ;
	     }
	     else{
	        die if(/>/);
		    ## remove spaces 
		    s/ //g;
		    chomp;
		    $sequence = $sequence . $_ ;
	     }
}

## now print the sequences
PrintSequence($NAME,$sequence,$docheck);
close($ifh);

## have a first parse that save the length...
sub ProcessSingleSequence{
	 my ($name,$sequence,$docheck) = @_ ;
	 my $len = length($sequence) ;

	 my @l = split "_" , $name ;
	 my $N = @l -1 ;
	 $info->{$l[0]} = {} if(! defined $info->{$l[0]} );
	 $info->{$l[0]}->{$l[$N]} = $len ;
}

sub PrintSequence{
	 my ($name,$sequence,$docheck) = @_ ;
	 my @l = split "_" , $name ;
	 my $N = @l -1 ;
	 my $ID = $l[0] ;
	 my $NUM = $l[$N] ;
	 die if(!exists $best->{$ID});
	 my $tab = $best->{$ID};
	 if(exists $tab->{$NUM}){
	 	 my $RANK = $tab->{$NUM}; 
	 	 my $NMM = "$ID.ORF_$NUM";
		 my $OFHLISTSINGLE = util_open_or_append("FASTADIR_ORF/$ID.list");
		 print $OFHLISTSINGLE "$NMM\n";
         my $ofh = util_write("FASTADIR_ORF/$NMM.ALL.1.fasta");
	 	 print $ofh ">$NMM\n";
	 	 print $ofh "$sequence\n";
		 my $len = length($sequence);
		 close($ofh) ;
		 print $ofhlist "$NMM\n";

		 if($RANK eq 1){
	 	      print $ofhlistlongest ">$ID $NMM\n";
	 	      print $ofhlistlongest "$sequence\n";
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
