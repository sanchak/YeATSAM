
package AAConfig;
use MyUtils;
use PDB;
use Atom;
use Residue;
use MyGeom;
require Exporter;
@ISA = qw(Exporter );
#@EXPORT = qw($fields);

use strict ;
use Carp ;
use FileHandle ;
use Getopt::Long;
use vars qw($AUTOLOAD);
use Getopt::Long;
use File::Basename;

my $MINDIST = 1 ; 
my $statesseen = {};

my $fields = {
    POLAR => undef, 
    THREE2ONE => undef, 
    INDICES => undef, 
    SIMNODES => undef, 
    PATH => undef, 
    ATOMIDX => undef, 
    CONFIG => undef, 
};

sub new{
    my $that = shift ; 
    my $class = ref($that) || $that ;
	my ($infile) = @_ ;
    my $ifh = util_read($infile);
    my $self = {} ; 
    map { $self->{$_} = undef ; } (keys %{$fields});

	$self->{THREE2ONE} = {};
	$self->{POLAR} = {};

	while(<$ifh>){
	    next if(/^\s*$/);
	    next if(/^\s*#/);
		my ($nm,$three,$single,$polar,$basic,$acidic,$amide) = split ;
		$three = uc($three);
		$self->{THREE2ONE}->{$three} = $single ; 

		$self->{POLAR}->{$single} = $polar ; 
		$self->{BASIC}->{$single} = $basic ; 
		$self->{ACIDIC}->{$single} = $acidic ; 
		$self->{AMIDE}->{$single} = $amide ; 
	}

    bless $self, $class ; 
    $self ;
}

sub AUTOLOAD {
    my $self = shift;
    my $attr = $AUTOLOAD;
    $attr =~ s/.*:://;
    return unless $attr =~ /[^A-Z]/;  # skip DESTROY and all-cap methods
    croak "invalid attribute method: ->$attr()" unless exists $fields->{$attr} ; 
    $self->{$attr} = shift if @_;
    return $self->{$attr};
}

sub MakeSingle{
	my ($self,$nm) = @_ ; 
	croak "NM not defined" if(!defined $nm);
	if(length($nm) == 3){
		$nm = $self->{THREE2ONE}->{$nm} ;
	   croak "NM not defined after mapping" if(!defined $nm);
	}
	return $nm ; 
}

sub IsPolar{
	my ($self,$nm) = @_ ; 
	$nm = $self->MakeSingle($nm);
	return $self->{POLAR}->{$nm}; 
}
sub IsBasic{
	my ($self,$nm) = @_ ; 
	$nm = $self->MakeSingle($nm);
	return $self->{BASIC}->{$nm}; 
}
sub IsAcidic{
	my ($self,$nm) = @_ ; 
	$nm = $self->MakeSingle($nm);
	return $self->{ACIDIC}->{$nm}; 
}
sub IsAmide{
	my ($self,$nm) = @_ ; 
	$nm = $self->MakeSingle($nm);
	return $self->{AMIDE}->{$nm}; 
}
sub IsElectrostatic{
	my ($self,$nm1,$nm2) = @_ ; 
	$nm1 = $self->MakeSingle($nm1);
	$nm2 = $self->MakeSingle($nm2);
	my $ret = 0;
	if($self->IsBasic($nm1) && $self->IsAcidic($nm2)){
		$ret = 1 ;
	}
	elsif ($self->IsBasic($nm2) && $self->IsAcidic($nm1)){
		$ret = -1 ;
	}
	return $ret ;
}
