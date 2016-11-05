package KMER;
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
    NAME => undef, 
    SEQUENCE => undef, 
    KSIZE => undef, 
    SEQUENCEFILE => undef ,
    ISNT => undef ,
    MASKL => undef ,
    SLIDINGUNIT => undef ,
    LEN => undef ,
};

sub new{
    my $that = shift ; 
    my $class = ref($that) || $that ;

    my $self =  {};

    map { $self->{$_} = undef ; } (keys %{$fields});

	$self->{NAME} = shift ;
	$self->{SEQUENCE} = shift ;
	$self->{KSIZE} = shift ;
	$self->{isNT} = shift ;

    $self->{LEN} = length($self->{SEQUENCE});

	$self->{SLIDINGUNIT} = 1 ;

	$self->{SPLITSEQ} = {};


    bless $self, $class ; 
	#$self->AddResidue($self->{PSEUDORESIDUE}); 
    $self ;
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
sub SetMASKL{
    my ($self,$what) = @_ ;
	## HACK for plants
	$self->{SEQUENCE} =~ s/$what//g ;
    $self->{LEN} = length($self->{SEQUENCE});
	if($self->{LEN} < $self->{KSIZE}){
		return -1 ;
	}
}

sub SanityCheck{
    my ($self) = @_ ;

	if($self->{LEN} < $self->{KSIZE}){
		return -1 ;
	}
	return 0 ;
}


sub SplitSliding{
    my ($self) = @_ ;

    return util_SplitSliding ($self->{SEQUENCE} , $self->{KSIZE}, $self->{NAME}, $self->{SLIDINGUNIT},$self->{SPLITSEQ} );
	
}


sub Print{
    my ($self) = @_ ;
	print "Original = $self->{SEQUENCE}\n";
	foreach my $key  (sort keys %{$self->{SPLITSEQ}}){
		my $v = $self->{SPLITSEQ}->{$key};
		print "$key $v \n";
	}
}


