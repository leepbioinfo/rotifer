# POD documentation - main docs before the code

=head1 NAME

Rotifer::SequenceModel - holder of annotation for a sequence 

=head1 SYNOPSIS

  # Simple example for single value (scalar) option
  use Rotifer::SequenceModel;
  my $profile = Rotifer::SequenceModel->new(name => "2C");
  print $profile->id;

=head1 DESCRIPTION

This module represents individual profiles in a Rotifer compatible
profile database and provides some data and facilities for
handling sequence profiles.

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 DEPENDENCIES

=over

=item Moose

=item File::Basename

=back

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Rotifer::SequenceModel;

# Load dependencies and pragmas
use Moose;
use File::Basename;
use Digest::MD5;
use strict;
use warnings;

our $VERSION = "0.1";

=head2 ATTRIBUTES

This Moose based class has the following attributes, all of which can
be set by passing the "attribute => value" to the new() or using the
attribute's accessor, which is a method with the same name, e.g.
$instance->attribute('new value').

=over

=item * name : this attribute is obligatory and corresponds to the user
defined name of a profile

=item * id   : the main identifier of each profile. This attribute 
should be unique, i.e. should not be used more than once in any 
analysis.

=item * file : fully qualified path to the profile itself or to a
database that contains the profile

=item * extension : file name extension

=item * number_of_iterations : number of rounds in iterated searches

=item * inclusion_threshold  : maximum e-value for inclusion in the
HMM or PSSSM in the next round of an iterated search

=item * terms => hash of associated annotations (term => type)

=back

=cut

# Name is obligatory
has 'name' => (is       => 'rw',
	       isa      => 'Str',
	       required => 1,
	       trigger  => sub {
		   my ($self, $new, $old) = @_;
		   die "ERROR: profile names cannot include any of the following characters:
 backslash (\\) comma (,) pipe (|) semicolon (;) slash (/) interrogation (?)
and may contain only one dot (.) near the end of the name. 
If present, this dot must be followed by a number that indicates
the use of a specific seed query or alignment to build the profile.
Conflicting profile: $new" if ($new =~ /([\,\/\|\\\;\?])(\.\d*)?$/);
	       },
    );

# ID is uc(NAME) unless defined by the user
has 'id' => (is       => 'rw',
	     isa      => 'Str',
	     lazy     => 1, # Wait until name is defined
	     default  => sub { my $self = shift; return uc($self->name) }
    );

# DOMAIN: same as name, unless specified
has 'annotation' => (is        => 'rw',
		     isa       => 'HashRef[ArrayRef[ArrayRef|Object|Str]]',
		     traits    => [ 'Hash' ],
		     handles   => {
			 has_annotation           => 'exists',
			 add_Annotation           => 'set',
			 get_Annotations          => 'get',
			 get_all_annotation_keys  => 'keys',
			 remove_annotation        => 'delete',
		     },
		     default   => sub { {} }
    );

# File that contains the profile
has 'file' => (is  => 'rw',
	       isa => 'Str|Undef',
	       default => undef,
	       trigger => sub {
		   my ($self, $new, $old) = @_;
		   return 1 unless (defined $new);
		   my $ext = $self->extension; 
		   die "Unable to read profile ".$self->name." ($new)" unless (-r $new);
           #if ($ext =~ /^\.(chk|asn1)$/) {
           #    (my $query = $new) =~ s/\.(chk|asn1)$/.csq/;
           #    die "Cannot load query sequence file $query for profile $new" if ( ! -e $query );
           #}
		   return 1;
	       },
    );

# File checksum
has 'checksum' => (is  => 'rw',
		   isa => 'Str|Undef',
		   lazy => 1,
		   default => undef,
		   trigger => sub {
		       my ($self, $new, $old) = @_;
		       return 1 unless (defined $new);
		       my $ext = $self->file;
		       die "Unable to read profile ".$self->name." ($new)" unless (-r $new);
		       if ($ext =~ /^\.(chk|asn1)$/) {
			   (my $query = $new) =~ s/\.(chk|asn1)$/.csq/;
			   die "Cannot load query sequence file $query for profile $new" if ( ! -e $query );
		       }
		       return 1;
		   },
    );

# Run parameters
has 'executable'            => (is => 'rw', isa => 'Str', default => '');
has 'fetcher'               => (is => 'rw', isa => 'Str', default => '');
#has 'number_of_processors'  => (is => 'ro', isa => 'Int', default => \&_detect_processors);
has 'number_of_iterations'  => (is => 'rw', isa => 'Int', default => 3);
has 'inclusion_threshold'   => (is => 'rw', isa => 'Str', default => 0.01);
has 'maximum_evalue_cutoff' => (is => 'rw', isa => 'Str', default => 0.01);

=head2 PROFILE ANNOTATION

=head2 domain

 Title   : domain
 Usage   : $profile->domain()
 Function: retrieve the domain name of this profile
 Returns : string or undef
 Args    : none

=cut

sub domain {
    my ($self, $name) = @_;
    if (defined $name) {
	$self->remove_annotation('domain') if ($self->has_annotation('domain')); # Enforce one domain per profile
	$self->add_terms('domain',$name);
    }
    return undef unless $self->has_annotation('domain');
    my @domain = $self->get_terms('domain');
    return wantarray ? @domain : $domain[0];
}

=head2 add_terms

 Title   : add_terms
 Usage   : $profile->add_terms($category, @terms)
 Function: associate annotation with this profile
 Returns : number of terms added
 Args    : annotation category and terms

=cut

sub add_terms {
    my ($self, $category, @terms) = @_;
    if ($self->has_annotation($category)) {
	push(@{ $self->annotation->{$category} }, @terms);
    } else {
	$self->add_Annotation($category => [ @terms ]);
    }
    return scalar(@terms);
}

=head2 get_terms

 Title   : get_terms
 Usage   : $profile->get_terms($category)
 Function: retrieve all annotations terms
 Returns : array of strings (maybe empty)
 Args    : (optional) category

NOTE: when called without arguments, this method will
      retrieve all annotation terms but the user won't
      to which category each term belongs.

=cut

sub get_terms {
    my ($self, $category) = @_;
    my @keys = $self->get_all_annotation_keys;
    if (defined $category) {
	return () unless ($self->has_annotation($category));
	@keys = ($category);
    }
    return map { map { ref($_) eq "ARRAY" ? @{$_} : $_ } @{ $self->get_Annotations($_) } } @keys;
}

=head2 remove_terms

 Title   : remove_terms
 Usage   : $profile->remove_terms($category)
 Function: remove annotation terms
 Returns : nothing
 Args    : (option) annotation category

This method will remove all annotation if a category
is not specified

=cut

sub remove_terms {
    my ($self, $category) = @_;
    if (defined $category) {
	$self->remove_annotation($category);
    } else {
	@{ $self->annotation } = {};
    }
    return 1;
}

=head2 FILE HANDLING

=head2 extension

 Title   : extension
 Usage   : $profile->extension()
 Function: retrieve profile filename extension, e.g. .hmm
 Returns : string or undef
 Args    : none

=cut

sub extension {
    my $self = shift;
    return undef unless (defined $self->file);
    $self->file =~ /(\.[^\.]+)$/;
    return $1 || '';
}

=head2 fetch

 Title   : fetch
 Usage   : $profile->fetch()
 Function: retrieve a text representation of the profile
 Returns : string or undef (if cannot dump the profile)
 Args    : (optional) library path

This methods returns 

=cut

sub fetch {
    my ($self) = @_;

    # Try library when file is not set
    my $file = $self->file;
    return undef if (!defined $file);
    return undef if (!defined $self->extension);

    # Select retrieval method
    my $string = undef;
  METHOD: {
      my $ext = $self->extension;
      $ext =~ /\.chk|\.asn1/ && do {
	  (my $csq = $file) =~ s/\.chk$/.csq/;
	  open(PF,"<$csq") || die " Could not read profile $file";
	  $string = join("",<PF>);
	  close(PF);
	  last METHOD;
      };

      $ext =~ /\.h3i|\.ssi/ && do {
	  $file =~ s/${ext}$//;
	  my $name = $self->id;
	  my $prog = $self->fetcher;
	  return undef unless (defined $prog);
	  open(PF,"$prog $file $name 2>> /dev/null |") || die " Could not read library $file";
	  $string = join("",<PF>);
	  close(PF);
	  last METHOD;
      };

      $ext =~ /\.rps|\.pal/ && do {
	  my ($filename, $directories, $suffix) = fileparse($self->file, '.rps', '.pal');
	  return "RPS-BLAST $filename";
      };

      $ext =~ /\.fa|\.seq|\.fasta|\.hmm/ && do {
	  open(PF,"<$file") || die " Could not read profile $file";
	  $string = join("",<PF>);
	  close(PF);
	  last METHOD;
      }
    }

    return length($string) ? $string : undef;
}

# Moose-style library ending (to make things work and be fast and clean)
no Moose;
__PACKAGE__->meta->make_immutable;
