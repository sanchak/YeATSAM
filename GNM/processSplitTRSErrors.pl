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
my ($fastafile,$idx,$infile,$dirout,$dirin,$postfix,$p1,$p2,$outfile,$cutoff,$which_tech,$listfile,$protein);
my ($ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 1 ;
GetOptions(
            "dirout=s"=>\$dirout,
            "dirin=s"=>\$dirin ,
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
usage( "DOnt give option -outfile ") if(defined $outfile);
$outfile = "$infile.parsed";
my $ofhsame = util_write("$outfile.same");
my $ofhdiff = util_write("$outfile.diff");
my $ofhonlyone = util_write("$outfile.onlyone");
my $ofhexacttrs = util_write("$outfile.annotated");
my $ofhmapnames = util_write("$outfile.mapnames");
my $ofhcommands = util_write("$outfile.commands.csh");
my $ofhlist = util_write("$outfile.list");
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a input file name => option -dirin ") if(!defined $dirin);
usage( "Need to give a input file name => option -dirout ") if(!defined $dirout);
my $ifh = util_read($infile);

my $fastaext = "ALL.1.fasta";

my $info = {};
my $fullinfo = {};
my $savePerTRS = {};
while(<$ifh>){
     next if(/^\s*$/);
	 next if(/^\s*#/);

     chomp ;
	 my ($nm,$junk) = split ; 
	 my $orig = $nm ;
	 $junk =~ s/\..//;
	 $info->{$nm} = $junk ;
	 $fullinfo->{$nm} = $_ ;

	 $nm =~ s/.ORF.*//;
	 $savePerTRS->{$nm} = [] if(!defined $savePerTRS->{$nm});
	 push @{$savePerTRS->{$nm}}, $orig;
}

foreach my $nm (sort keys %{$savePerTRS}){
	print $ofhexacttrs "$nm\n";
	my @lll = @{$savePerTRS->{$nm}};
	my $tabletmp = {};
	foreach my $i (@lll){
	    my $a = $info->{$i};
		$tabletmp->{$a} = [] if(! defined $tabletmp->{$a});
		push @{$tabletmp->{$a}}, $i ;
	}
	my @KEYS = (keys %{$tabletmp});
	my $NNNN = @KEYS;
	my @l;
	if($NNNN eq 1){
		my $a = $KEYS[0];
		@l = @{$tabletmp->{$a}};
	}
	else{
		foreach my $a (@KEYS){
			my @q = @{$tabletmp->{$a}};
			my $first = $q[0];
			push @l, $first ;
		}
	}


    @l = sort @l ;
	if(@l eq 1){
		my $orf = $l[0];
		my $a = $fullinfo->{$l[0]};
		print $ofhonlyone "$nm @l $a\n";
		print $ofhlist "$nm\n";
		print $ofhcommands "\\cp -f $dirin/$orf.$fastaext $dirout/$nm.$fastaext #uniq\n";
		print $ofhmapnames "$nm $orf\n";
		next ;
	}

	my $fullanno = $fullinfo->{$l[0]};
	my $a = $info->{$l[0]};
	my $b = $info->{$l[1]};
	if(@l eq 2){
			 my $orf1 = $l[0];
			 my $orf2 = $l[1];
			 my $F1 = "$dirin/$orf1.$fastaext";
			 my $F2 = "$dirin/$orf2.$fastaext";
			 my ($str1) = util_readfasta($F1);
			 my ($str2) = util_readfasta($F2);
			 my $LEN1 = length($str1);
			 my $LEN2 = length($str2);
		if($a eq $b){
			 my $FINALORF = $LEN1 > $LEN2 ? $F1 : $F2 ;
			 my $orfname = $LEN1 > $LEN2 ? $orf1 : $orf2 ;
		     print $ofhcommands "\\cp -f $FINALORF $dirout/$nm.$fastaext #same2 choosing between $LEN1 and $LEN2\n";
		     print $ofhmapnames "$nm $orfname\n";
		     print $ofhlist "$nm\n";

		     print $ofhsame "$nm 2 $fullanno @l $a $b \n";
		}
	    else{
			 my $nmA = "$nm" . "_A";
			 my $nmB = "$nm" . "_B";
		     print $ofhcommands "\\cp -f $F1 $dirout/$nmA.$fastaext #diff $LEN1 and $LEN2\n";
		     print $ofhcommands "\\cp -f $F2 $dirout/$nmB.$fastaext #diff $LEN1 and $LEN2\n";
		     print $ofhmapnames "$nmA $orf1 \n";
		     print $ofhmapnames "$nmB $orf2 \n";
		     print $ofhlist "$nmA\n";
		     print $ofhlist "$nmB\n";

		     print $ofhdiff "$nm 2 $fullanno @l $a $b \n";
		}
	}
	else{
			 my $orf1 = $l[0];
			 my $orf2 = $l[1];
			 my $orf3 = $l[2];
			 my $F1 = "$dirin/$orf1.$fastaext";
			 my $F2 = "$dirin/$orf2.$fastaext";
			 my $F3 = "$dirin/$orf3.$fastaext";
			 my ($str1) = util_readfasta($F1);
			 my ($str2) = util_readfasta($F2);
			 my ($str3) = util_readfasta($F3);
			 my $LEN1 = length($str1);
			 my $LEN2 = length($str2);
			 my $LEN3 = length($str3);
		my $c = $info->{$l[2]};
		if($a ne $b &&  $b ne $c && $a ne $c){ ## all three diff
			 my $nmA = "$nm" . "_A";
			 my $nmB = "$nm" . "_B";
			 my $nmC = "$nm" . "_C";
		     print $ofhcommands "\\cp -f $F1 $dirout/$nmA.$fastaext #diff3 $LEN1 \n";
		     print $ofhcommands "\\cp -f $F2 $dirout/$nmB.$fastaext #diff3 $LEN2 \n";
		     print $ofhcommands "\\cp -f $F3 $dirout/$nmC.$fastaext #diff3 $LEN3 \n";
		     print $ofhmapnames "$nmA $orf1 \n";
		     print $ofhmapnames "$nmB $orf2 \n";
		     print $ofhmapnames "$nmC $orf3 \n";
		     print $ofhlist "$nmA\n";
		     print $ofhlist "$nmB\n";
		     print $ofhlist "$nmC\n";

		     print $ofhdiff "$nm 3 $fullanno @l $a $b $c \n";
		}
		elsif($a eq $b &&  $b eq $c){ ## all three same 
			 my $FINALORF = $LEN1 > $LEN2 ? $F1 : $F2 ;
			 my $ORFNAME = $LEN1 > $LEN2 ? $orf1 : $orf2 ;
			 my $FINALLEN = $LEN1 > $LEN2 ? $LEN1 : $LEN2 ;
			 $FINALORF = $FINALLEN > $LEN3 ? $FINALORF : $F3 ;
			 $ORFNAME = $FINALLEN > $LEN3 ? $ORFNAME : $orf3 ;

		     print $ofhcommands "\\cp -f $FINALORF $dirout/$nm.$fastaext #same3 choosing between $LEN1 and $LEN2\n";
			 print $ofhmapnames "$nm $ORFNAME\n";
		     print $ofhsame "$nm 3 $fullanno @l $a $b $c \n";
		     print $ofhlist "$nm\n";
		}
		else{
			 die ;
		}
	}
}

system ("wc -l $outfile.*");

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
