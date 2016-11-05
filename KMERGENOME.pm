package KMERGENOME ;
use MyUtils;
use PDB;
use ConfigPDB;
use Atom;
use Residue;
use MyGeom ;
use Algorithm::Combinatorics qw(combinations) ;
use Math::Geometry ; 
#use Math::Trig;
#use Math::Trig ':radial';
require Exporter;
@ISA = qw(Exporter );
@EXPORT = qw($fields);
use Math::VectorReal qw(:all);  # Include O X Y Z axis constant vectors

use strict ;
use Carp ;
use FileHandle ;
use Getopt::Long;
use vars qw($AUTOLOAD);
use Getopt::Long;
use Cwd ;
use File::Basename;

no warnings 'redefine';


my $verbose = 0 ;
my $debugone = 0 ;

my $fields = {
    NAMES => undef, 
    SEQUENCES => undef, 
    KSIZE => undef, 
    SEQUENCEFILE => undef ,
    SPLITSEQ => undef ,
    ISNT => undef ,
    MASKL => undef ,
    SLIDINGUNIT => undef ,
    HOPPING => undef ,
    REVERSE => undef ,
    DELFACTOR => undef ,
    LEN => undef ,
    SAVEINDIVIDUAL => undef ,
};

sub new{
    my $that = shift ; 
    my $class = ref($that) || $that ;

    my $self =  {};

    map { $self->{$_} = undef ; } (keys %{$fields});

	#$self->{NAME} = shift ;
	$self->{SEQUENCEFILE} = shift ;
	$self->{KSIZE} = shift ;
	$self->{HOPPING} = shift ;
	$self->{isNT} = shift ;
	$self->{SAVEINDIVIDUAL} = shift ;
	$self->{REVERSE} = shift ;
	$self->{DELFACTOR} = shift ;


	$self->{SEQUENCES} = {};




	$self->{SLIDINGUNIT} = 1 ;
	$self->{SPLITSEQ} = {};
	$self->{NAMES} = {};


    bless $self, $class ; 

	$self->Split();
	
    $self ;
}


sub Split{
	my ($self) = @_ ;
	if($self->{HOPPING}){
		 print "Genome: hopping - simple !!\n";
    }
	else{
		 print "Genome: Sliding - careful !!\n";
	}
    my $ifh = util_read($self->{SEQUENCEFILE});
	my $NAME ;
	my $STR ;
    while(<$ifh>){
	     next if(/^\s*$/);
         if(/^\s*>/){
	 	    ## process one fasta
	 	    if(defined $NAME){
	           #$self->_ProcessSingleFasta($NAME,$STR,$expr2ignore);
			   ## expr2ignore is too long for AA 
	           $self->_ProcessSingleFasta($NAME,$STR);
               $STR = "";
		    }
    
		    s/>//;
		    ($NAME) = split ;
			die "Error: repeat name $NAME" if(exists $self->{NAMES}->{$NAME});
			$self->{NAMES}->{$NAME} = 1 ;
	     }
	     else{
	        die if(/>/);
		    ## remove spaces 
		    s/ //g;
		    chomp;
		    $STR = $STR . $_ ;
	     }
    }
	#$self->_ProcessSingleFasta($NAME,$STR,$expr2ignore);
	$self->_ProcessSingleFasta($NAME,$STR);
    close($ifh);
}


my $docheck = 1 ;
sub _ProcessSingleFasta{
    my ($self,$name,$sequence) = @_ ;
	my $len = length($sequence);
	if($docheck && $len > 60 ){
		$docheck = 0 ;
		util_CheckSequence($name,$sequence,$self->{isNT});
	}

	if($self->{REVERSE}){
	     $sequence = util_getComplimentaryString($sequence) ;
	}


	die if(!defined $name );
	if($self->{HOPPING}){
         util_SplitHopping ($sequence, $self->{KSIZE},$self->{DELFACTOR}, $name,$self->{SPLITSEQ});
	}
	else{
         util_SplitSliding ($sequence, $self->{KSIZE},$name, $self->{SLIDINGUNIT},$self->{SPLITSEQ});
	}

	if($self->{SAVEINDIVIDUAL}){
		$self->{SEQUENCES}->{$name} = $sequence ;
	}

}

sub AUTOLOAD {
  my $self = shift;
   my $type = ref ($self) || croak "$self is not an object";
   my $field = $AUTOLOAD;
   $field =~ s/.*://;
   unless (exists $self->{$field})
   {
      croak "$field does not exist in object/class $type";
   }
   if (@_)
   {
      $self->{$field} = shift;
   }
   else
   {
      return $self->{$field};
   }
}



sub FREE{
    my ($self) = @_ ;
	print "Freeing genome hash\n";
	foreach my $key  (sort keys %{$self->{SPLITSEQ}}){
		    undef $self->{SPLITSEQ}->{$key};
	}
	undef %{$self->{SPLITSEQ}};
}

sub Print{
    my ($self,$verbose) = @_ ;
	if($verbose){
	    foreach my $key  (sort keys %{$self->{SPLITSEQ}}){
		    my $v = $self->{SPLITSEQ}->{$key};
		    print "$key $v \n";
	    }
	}
	my $Nsequences = keys %{$self->{NAMES}};
	my $Nsplitsequences = keys %{$self->{SPLITSEQ}};
	print "GenomeInfo: There are $Nsequences sequences which are split in $Nsplitsequences with kmer size $self->{KSIZE}\n";
}


sub GetAllSequences{
    my ($self) = @_ ;
	return $self->{SEQUENCES};
}
