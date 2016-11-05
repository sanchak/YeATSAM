#!/usr/bin/perl -w 
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($all,$infile,$outfile,$or,$silent,$groupinfo);
my ($fastadir,$DIR,$listfile,$ignorefile,$mergedir);
my $size ;
my $cutofflength ;
my @types = (); 
my $annofile ;
my @motifs = (); 
GetOptions(
            "all"=>\$all ,
            "groupinfo"=>\$groupinfo ,
            "silent"=>\$silent ,
            "infile=s"=>\$infile ,
            "fastadir=s"=>\$fastadir ,
            "dir=s"=>\$DIR ,
            "mergedir=s"=>\$mergedir ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "size=i"=>\$size ,
            "or=i"=>\$or ,
            "cutofflength=i"=>\$cutofflength ,
            "type=s"=>\@types,
            "annofile=s"=>\$annofile,
            "motif=s"=>\@motifs,
            "outfile=s"=>\$outfile 
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);
usage( "Need to give a input file name => option -mergedir ") if(!defined $mergedir);
usage( "Need to give a protein -option -annofile  ") if(!defined $annofile);
usage( "Need to give a fastadir -option -fastadir  ") if(!defined $fastadir);
usage( "Need to give a protein -option -size  ") if(!defined $size);
$infile = $infile. $size ;
$outfile = "$infile.out";
my $ofh = util_write("$mergedir/$outfile");
my $ofhmapping = util_write("$mergedir/$outfile.mapping");
my $ofhlist = util_open_or_append("merge.list");

my $ofhdonealready = util_open_or_append("alreadydone");

my $alreadydone = {};
my @lll= util_read_list_sentences("alreadydone");
map { s/\s*//g ; $alreadydone->{$_} = 1 ; } @lll ;


my $ignoretable = {};
if(defined $ignorefile){
   my @lll= util_read_list_sentences($ignorefile);
   map { s/\s*//g ; $ignoretable->{$_} = 1 ; } @lll ;
}


my $ifh = util_read("$mergedir/$infile");
my $ENDDONE = {};
my $BEGINDONE = {};
my $STRINGS = {};

my $notprinted = 1 ;
my $seen = {};
while(<$ifh>){
	my @l = split ;
	if(@l eq 3){
		my ($j,$begin,$end) = @l ;
		die "Dont expect this because of the previous step" if(exists $BEGINDONE->{$begin});
                $BEGINDONE->{$begin} = 1 ;

		my @l = split "",$j ;
		my $letters = {};
		foreach my $l (@l){
			$letters->{$l} = 1 ;
		}
		my $N = (keys %{$letters});
		if($N < 3){
			if(!exists $seen->{$j}){

			my $ofhignoremerged = util_open_or_append("ignoremerged");
			print $ofhignoremerged "$end $begin xx_rep=$j\n";

			   print "Seems like repetitive: " if($notprinted) ;
			   $seen->{$j} = 1 ;
			   $notprinted = 0 ;
			   print "-- $j -- ";
			}
			next ;
		}
		else{
		    $ENDDONE->{$end} = [] if(!defined $ENDDONE->{$end});
		    push @{$ENDDONE->{$end}}, $begin;
			$STRINGS->{$end} = $j ;
		}
	}
}
print "\n";


my $ifhann = util_read($annofile);
my $infoann = {};
while(<$ifhann>){
    next if(/Warning/i);
     my ($nm,$junk) = split ;
     $nm =~ s/.blast.nt//;	
     $infoann->{$nm} = $junk ;
}



my $ofhwarning = util_write("$mergedir/warning.overlap.$size");

my $ofhinclude = util_write("$mergedir/list.$infile.include");
my $ofhexclude = util_write("$mergedir/list.$infile.exclude");

### Choose only those which have one to one 
my $CCCDiff  = 0 ;
my $CCCSame  = 0 ;
foreach my $end (keys %{$ENDDONE}){
	my @l = @{$ENDDONE->{$end}};
	my $SAME = 0 ;
	if(@l eq 1){
		my $begin = $l[0];

	   my $str = $STRINGS->{$end};

		my $JOINEDSTR = $begin.$end ;
		if($JOINEDSTR =~ /_(A|B|C)/){
			my $ofhignoremerged = util_open_or_append("ignoremerged");
			print $ofhignoremerged "$end $begin xx_merged=$str\n";
		}




		my ($scaff1) = ($begin =~ /C(.*)\_G/);
		my ($scaff2) = ($end =~ /C(.*)\_G/);

		my $SCAFFDIFF = 0 ;
		if($scaff2 ne $scaff1){
		    print $ofh "DIFF $end $begin $str\n";
	            $CCCDiff++;
		    $SCAFFDIFF = 1 ;
		}

		{
		  if(0){
                   $end =~ s/_A/_a/;
                   $end =~ s/_B/_b/;
                   $end =~ s/_C/_c/;

                   $begin =~ s/_A/_a/;
                   $begin =~ s/_B/_b/;
                   $begin =~ s/_C/_c/;
		   }

		   if(exists $alreadydone->{$end} || exists $alreadydone->{$begin}){
		   	next ;
                   }
		   print $ofhdonealready "$end\n";
		   print $ofhdonealready "$begin\n";
		   $alreadydone->{$end} = 1;
		   $alreadydone->{$begin} = 1;

			my $A = $end ;
			my $B = $begin ;
			$A =~ s/.ORF.*//;
			$B =~ s/.ORF.*//;

			my $GENOMEA = $infoann->{$A} or die " $A ";
			my $GENOMEB = $infoann->{$B} or die "$B";
			my $DIFFGENOME = 0 ;
			if($GENOMEB ne $GENOMEA){
				print $ofhwarning "WARN $end has $GENOMEA, while $begin has $GENOMEB\n";
		                print $ofh "DIFFGENOME $end $begin $str\n";
				$DIFFGENOME = 1;
			}

			if($DIFFGENOME && $SCAFFDIFF && $size < 6 ){
			}
			else{
					my $newnm = "$end" . "_MERGED";
		            print $ofh "SAME $end $begin $str\n";
		            print $ofhlist "$end $begin $newnm\n";
			    $CCCSame++;
			    $SAME = 1 ;
			}
		}

		# We write for all now
		if(1){
           $end =~ s/;//g;
           $begin =~ s/;//g;
		   #my $newnm = "$end.$begin";
		   my $newnm = "$end" . "_MERGED";

		   if(0){
                   $end =~ s/_A/_a/;
                   $end =~ s/_B/_b/;
                   $end =~ s/_C/_c/;

                   $begin =~ s/_A/_a/;
                   $begin =~ s/_B/_b/;
                   $begin =~ s/_C/_c/;
		   }
		   my ($fe,$x) =  util_readfasta("$fastadir/$end.ALL.1.fasta");
		   my ($fb,$y) =  util_readfasta("$fastadir/$begin.ALL.1.fasta");

		   foreach my $i (1..$size){
		       $fb =~ s/.//;
		   }

		   my $FH = util_write("$fastadir/$newnm.ALL.1.fasta");
		   my $join = $fe . $fb ;
		   print $FH ">$newnm\n";
		   print $FH "$join\n";

		   if(0) # for nucleotide
		   {
               $end =~ s/;//g;
               $begin =~ s/;//g;
			   $end =~ s/.OR.*//;
			   $begin =~ s/.OR.*//;
		   my ($fe,$x) =  util_readfasta("FASTADIR_NT/$end.ALL.1.fasta");
		   my ($fb,$y) =  util_readfasta("FASTADIR_NT/$begin.ALL.1.fasta");

		   my $FH = util_write("FASTADIR_NT/$newnm.ALL.1.fasta");
		   my $join = $fe . $fb ;
		   print $FH ">$newnm\n";
		   print $FH "$join\n";
		   close ($FH);
		   }


		   if($SAME){
		       print $ofhinclude "$newnm\n";


                       $end =~ s/;//g;
                       $begin =~ s/;//g;

		       print $ofhexclude "$end\n";
		       print $ofhexclude "$begin\n";
		   }

		}
	}
}

print "found $CCCSame matches with same prefix and $CCCDiff matches with different prefixes  \n";

exit ;


sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
