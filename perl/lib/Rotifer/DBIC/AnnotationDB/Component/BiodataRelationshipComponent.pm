=head1 NAME

Rotifer::DBIC::AnnotationDB::Component::BiodataRelationshipComponent - methods for biosequences

=head1 DESCRIPTION

This modules extends Rotifer::DBIC::Result::BiodataRelationship.

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

package Rotifer::DBIC::AnnotationDB::Component::BiodataRelationshipComponent;

use strict;
use warnings;
use Carp::Clan qr/^Rotifer::DBIC::AnnotationDB/;
use Moose;
extends qw(
 Rotifer::DBIC::AnnotationDB::Component::WithDataset
 DBIx::Class
);

=head2 length

 Title   : length
 Usage   : $text = $rs->length()
 Function: get the sequence length
 Returns : integer or undef
 Args    : none

=cut

sub length {
    my $self = shift;
    return $self->end - $self->start_coord + 1 if (defined $self->end_coord && defined $self->start_coord);
    return undef;
}

__PACKAGE__->meta->make_immutable(inline_constructor => 0);
1;
