#!/usr/bin/perl -w 
use strict ;
use FileHandle ;
use Getopt::Long;
use Cwd ;
use MyUtils;
use MyGeom;
use PDB;
use ConfigPDB;
use Math::Geometry ;
use Math::Geometry::Planar;


use Time::HiRes qw( usleep ualarm gettimeofday tv_interval clock_gettime clock_getres  clock);
use POSIX qw(floor);
my $commandline = util_get_cmdline("",\@ARGV) ;
my ($finalvalue,$infile,$p1,$p2,$outfile,$cutoff,$which_tech,$map2scaffold,$listfile,$map2len);
my ($ignorefile,@expressions);
my $howmany = 100000 ;
my $verbose = 0 ;
GetOptions(
            "which_tech=s"=>\$which_tech ,
            "map2len=s"=>\$map2len ,
            "infile=s"=>\$infile ,
            "map2scaffold=s"=>\$map2scaffold ,
            "p2=s"=>\$p2 ,
            "listfile=s"=>\$listfile ,
            "ignorefile=s"=>\$ignorefile ,
            "outfile=s"=>\$outfile ,
            "expr=s"=>\@expressions,
            "howmany=i"=>\$howmany ,
            "finalvalue=i"=>\$finalvalue ,
            "verbose=i"=>\$verbose ,
            "cutoff=f"=>\$cutoff ,
           );
die "Dont recognize command line arg @ARGV " if(@ARGV);
usage( "Need to give a output file name => option -outfile ") if(!defined $outfile);
usage( "Need to give a input file name => option -infile ") if(!defined $infile);

usage( "Need to give a listfile -option -listfile  ") if(!defined $listfile);
usage( "Need to give a map2len pdb id -option -map2len  ") if(!defined $map2len);
usage( "Need to give a map2scaffold pdb id -option -map2scaffold  ") if(!defined $map2scaffold);
usage( "Need to give a ignorefile pdb id -option -ignorefile  ") if(!defined $ignorefile);
my $CNT = 0 ; 
#my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC) = util_SetEnvVars();
 my ($RESULTDIR,$PDBDIR,$FASTADIR,$APBSDIR,$FPOCKET,$SRC,$MATCH3D,$ANNDIR, $UNIPROT,$PREMONITION,$HELIXDIR,$DSSP,$CONFIGGRP,$BLASTOUT,$BLASTDB) = util_SetEnvVars();

my $PWD = cwd;

my ($map2lentable) = util_maketablefromfile($map2len);
my ($map2scaffoldtable) = util_maketablefromfile($map2scaffold);



my @list= util_read_list_words($listfile);
my $list = {};
map { s/\s*//g ; die if (exists $list->{$_}) ; $list->{$_} = 1 ; } @list ;



my @values = qw (2000 1000 500 250 125 70 );
my $NNN = @values -1 ;

if(!defined $finalvalue){
  $finalvalue = $values[$NNN];
}

my $OUTDIR = "INFO.$finalvalue";

$ignorefile = "$OUTDIR/$ignorefile";
unlink $ignorefile ;
system ("touch $ignorefile");

system ("mkdir -p $OUTDIR");
$outfile = "$OUTDIR/$outfile";
my $ofh = util_write($outfile);
my $SPLITLEVELSDIR = "$OUTDIR/SPLITLEVELS.$finalvalue";
system ("rm -rf $SPLITLEVELSDIR");
system ("mkdir -p $SPLITLEVELSDIR");
system ("mkdir -p $OUTDIR/CLUBBED");

my $info = {};
my $level = 0 ;
while(@values){
    my $ofhignore = util_open_or_append("$ignorefile");
    $level++;

    my $ITER = shift @values;
    last if(defined $finalvalue && $ITER < $finalvalue);
    
    my $comm = "groupBasedonCutoff.pl -outf $OUTDIR/CLUBBED/clubbed.$ITER -inf $infile -cutof $ITER";
    print "Info: Running $comm\n";
    system ($comm);

    my $ignoretable = {};
    my @lll= util_read_list_sentences($ignorefile);
    map { s/\s*//g ; $ignoretable->{$_} = 1 ; } @lll ;
    
    
    my $ifh = util_read("$OUTDIR/CLUBBED/clubbed.$ITER");
    while(<$ifh>){
         next if(/^\s*$/);
         next if(/^\s*#/);
    
         s/.*=N //;
         my (@l) = split ; 
         die if(@l < 2); ## cant be

         ## sort based on length of TRS
         my $sorted = {};
         foreach my $i (@l){
             next if(exists $ignoretable->{$i});
            die "$i does not exist in len" if (! exists $map2lentable->{$i});
             my $len = $map2lentable->{$i}; 
            #print "$i $len \n";
            $sorted->{$i} = $len ;

         }

         ## this can happen since we are ignoring...
         next if(keys %{$sorted} < 2);

         my @sl = sort { $sorted->{$b} <=> $sorted->{$a}} (keys %{$sorted});

         if(1){
         my $ADDEDONE = 0 ;
         foreach my $i (@sl){
            next if(!exists $info->{$i}) ;

             my $tmp = $sorted ;
            delete $tmp->{$i} ;
            my @ll = (sort keys %{$tmp});
            if(@ll){
               $ADDEDONE = 1;
               AddToOneLevel($ITER,$i,\@ll,$ignoretable,$ofhignore);
            }
         }
         next if($ADDEDONE);
         }

         my $cnt = 0 ;
         my $first ;
                 foreach my $i (@sl){
             ## Do not ignore first entry 
             if($cnt){
                if(exists $info->{$i}){
                    die "$first $i" ;
                }
             }
             else{
                 $first = $i;
                 #print "leng of $first = is $sorted->{$i} \n";
                              die if(exists $ignoretable->{$first});
             }
             $cnt++;
         }

         shift @sl ;
        
         if(! exists $info->{$first}){
             $info->{$first} = {};
         }
         else{
             ## just printing
            my @LL = (keys %{$info->{$first}});
            my $NN = @LL ;
            my $ADDTHISTIME = @sl ;
            print "INFO: $first was seen in previous iter with $NN values, adding $ADDTHISTIME \n" if($verbose);
        }

        AddToOneLevel($ITER,$first,\@sl,$ignoretable,$ofhignore);
    }
    close($ifh);
    close($ofhignore);

    ## info 
    my @TRS = (keys %{$info});
    my $NN = @TRS ;
    print "Found $NN in this iter $ITER \n";

}
sub AddToOneLevel{
    my ($ITER,$first,$SL,$ignoretable,$ofhignore) =@_ ;
    my @sl = @{$SL};
    my $fnm = "$SPLITLEVELSDIR/$first.level" ;
    my $prefix = (! -e $fnm) ? "$first\n" : "";
    my $ofhsplit = (! -e $fnm)? util_write($fnm) : util_append("$SPLITLEVELSDIR/$first.level");
    $prefix = $prefix . "\t $ITER ";
    print $ofhsplit "$prefix @sl \n";
    close($ofhsplit);

    foreach my $i (@sl){
           if(! exists $ignoretable->{$i}){
           print $ofhignore "$i\n";
       }
       $info->{$first}->{$i} = 1 ;
    }
}

my @TRS = (keys %{$info});
my $ofhscaffolds = util_write("$OUTDIR/mappedscaffolds");
my $scaffoldcnt = 0 ;
$, = " ";
print "=== Analyzing scaffolds ====\n";
foreach my $i (@TRS){

    my @l = (keys %{$info->{$i}});

    print $ofh "$i @l \n";

    ## include all TRS
    push @l, $i ;
    my $scaffolds = {};

    foreach my $x (@l){
        $x =~ s/_MERGED//;
        $x =~ s/_A//;
        $x =~ s/_B//;
        $x =~ s/_C//;
        my $scaffold = $map2scaffoldtable->{$x};

        die "Did not fined scaffold for $x" if(! defined $scaffold);
        $scaffolds->{$scaffold} = 1 ;
    }
    

    # Print all scaffolds mapped to the head TRS
    print $ofhscaffolds "$i ";
    foreach my $scaffold (keys %{$scaffolds}){
        print $ofhscaffolds " $scaffold ";
        $scaffoldcnt++;
    }
    print $ofhscaffolds " \n";
}
close($ofh);
#print "scaffoldcnt = $scaffoldcnt\n";

print "======= Unique gene list with at least one common gene ==== \n";
system("extractindexfromfile.pl -idx 0 -in $outfile");

my ($singlegeneWithCounterparts,$N1) = util_maketablefromfile($outfile);
my ($counterparts,$N2) = util_maketablefromfile($ignorefile);

print "singlegeneWithCounterparts = $N1, counterparts =$N2, scaffoldcnt = $scaffoldcnt \n";

#my ($mergetable) = util_table_merge($singlegeneWithCounterparts,$counterparts);
#my ($common,$inAbutnotinB,$inBbutnotinA) = util_table_diff($list,$mergetable);
#print "These are possible singletons $inAbutnotinB\n";

system ("mappingAddCount.pl -in $outfile -out $outfile.sort ");
system ("mappingAddCount.pl -in $OUTDIR/mappedscaffolds -out $OUTDIR/mappedscaffolds.sort ");
#system ("frequencyDistributionAbs.pl -outf freq -inf mappedscaffolds.sort -idx 0 -max 1000 -delta 1 -start 0");



my $tmpfile = "$outfile.mappedscaffolds";
unlink $tmpfile ;
print "cat $outfile.0 $ignorefile > \! $tmpfile\n";
print "TW $listfile $tmpfile\n";
print "extractthoseinlist.pl -in genome.annotated.csv -lis ofhinAbutnotinB -tag CCC\n";
print "extractthoseinlist.pl -in nonelist -lis ofhinAbutnotinB -tag CCC\n";

sub usage{
    my ($msg) = @_ ;
    print $msg , "\n" ; 
print << "ENDOFUSAGE" ; 
ENDOFUSAGE
    die ;
}
