=head1 NAME

Rotifer::DBIC::AnnotationDB::Component::BiosequenceComponent - methods for biosequences

=head1 DESCRIPTION
 
This modules extends Rotifer::DBIC::Result::Biosequence by overloading some
of DBIx::Class' default methods (update, insert, new) to automate handling of
the checksum column and to add support to the Bio::PrimarySeqI interface.

Checksum methods are expected to be replace by SQL triggers in the future.

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::DBIC::AnnotationDB::Component::BiosequenceComponent;

use strict;
use warnings;
use Carp::Clan qr/^Rotifer::DBIC::AnnotationDB/;
use Digest::SHA qw(sha256_hex);
use Moose;

#use base qw(DBIx::Class Bio::PrimarySeqI Bio::Root::Root);
extends qw(
 Rotifer::DBIC::AnnotationDB::Component::WithDataset
 DBIx::Class
 Moose::Object
 Bio::PrimarySeqI
 Bio::Root::Root
);

# Set default checksum to a SHA256 encoding of residues
sub new {
    my ( $class, $attrs ) = @_;

    # Create checksum
    $attrs->{checksum} = sha256_hex($attrs->{residues});
    $attrs->{seqlen}   = length($attrs->{residues}) || 0;

    return  $class->next::method($attrs);
}

# Set default checksum to a SHA256 encoding of residues
sub update {
    my ( $class, $attrs ) = @_;

    if (defined $attrs) {
	$attrs->{checksum} = sha256_hex($attrs->{residues})  if (! exists $attrs->{checksum} );
	$attrs->{seqlen}   = length($attrs->{residues}) || 0 if (! exists $attrs->{seqlen} );;
    }

    return  $class->next::method($attrs);
}

=head1 Bio::PrimarySeqI-compatible methods

=cut

=head2 display_id, primary_id, accession_number

 Mostly useless but implemented!

=cut

{ no warnings 'once';
  *display_id       = \&Rotifer::DBIC::AnnotationDB::Result::Biosequence::checksum;
  *primary_id       = \&Rotifer::DBIC::AnnotationDB::Result::Biosequence::biosequence_id;
  *accession        = \&Rotifer::DBIC::AnnotationDB::Result::Biosequence::checksum;
  *accession_number = \&Rotifer::DBIC::AnnotationDB::Result::Biosequence::checksum;
}

=head2 seq

 Title   : seq
 Usage   : $seq = $rs->seq(%attrs)
 Function: get/set residues for this sequence
 Returns : DBIx::Class::Row or undef
 Args    : hash reference

=cut

sub seq {
    return shift->residues(@_);
}

=head2 subseq

 Title   : subseq
 Usage   : $text = $rs->subseq(10,20)
 Function: Returns the subseq from start to end, where the first base
           is 1 and the number is inclusive, i.e. 1-2 are the first two
           bases of the sequence.

           Start cannot be larger than end but can be equal.

 Returns : string or undef
 Args    : pair of integers

=cut

sub subseq {
    my ($self, $start, $end) = @_;
    return undef if ($self->seqlen == 0);
    $start = 1 unless (defined $start);
    $end = $self->seqlen unless (defined $end);
    croak "fetching subsequences requires start <= end" if ($start > $end);
    my $length = $end - $start + 1;
    return substr($self->residues, $start-1, $length);
}

=head2 trunc

 Title   : trunc
 Usage   : $text = $rs->trunc(10,20)
 Function: same as subseq() but return a new object
 Returns : Rotifer::DBIC::AnnotationDB::Result::Biodata
 Args    : pair of integers

 Note    : only DBIx::Class::Row->new() is executed
           you must run $obj->insert() to store it

=cut

sub trunc {
    my ($self) = shift;
    return undef if ($self->seqlen == 0);
    my $seq = $self->subseq(@_);
    return $self->result_source->schema->resultset("Biosequence")->find_or_new({ term => $self->term, residues => $seq });
}

=head2 length

 Title   : length
 Usage   : $text = $rs->length()
 Function: get the sequence length
 Returns : integer or undef
 Args    : none

=cut

sub length {
    my $self = shift;
    return $self->seqlen || CORE::length($self->residues) || 0;
}

=head2 description

 Title   : description
 Usage   : $text = $rs->description()
 Function: does nothing
 Returns : undef
 Args    : 

=cut

sub description {
    return undef;
}

{ no warnings 'once';
  *desc = \&description;
};

=head2 alphabet

 Title   : alphabet
 Usage   : $text = $rs->alphabet()
 Function: does nothing
 Returns : undef
 Args    : 

=cut

sub alphabet {
    my $self = shift;
    my $term = $self->term->name;
    return $term eq "polypeptide" ? "protein" : $term;
}

=head2 is_circular

 Title   : is_circular
 Usage   : $text = $rs->is_circular()
 Function: find whether this sequence is a circular replicon
 Returns : boolean
 Args    : none

=cut

sub is_circular { shift->get_column('circular') };

# signal to BioPerl that this sequence can't be cloned
sub can_call_new { 0 }

# Make Perl and Moose happy and well behaved
__PACKAGE__->meta->make_immutable(inline_constructor => 0);
1;
