=head1 NAME

Rotifer::DBIC::AnnotationDB::Component::Annotatable - Bio::AnnotatableI compatibility

=head1 DESCRIPTION
 
This modules add the annotataion() methods to your class. You can use it to retrieve a
resultset that implements the Bio::Annotation::CollectionI interface.

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::DBIC::AnnotationDB::Component::Annotatable;

use strict;
use warnings;
use Carp::Clan qr/^Rotifer::DBIC::AnnotationDB/;
use base qw(DBIx::Class);

=head2 annotation

 Title   : annotation
 Usage   : $rs->annotation(%attrs)
 Function: get/set residues for this sequence
 Returns : DBIx::Class::Row or undef
 Args    : hash reference

=cut

sub annotation {
    my ($self) = shift;
    my $annotation_source = lc($self->result_source->source_name) . "_attribute_values";
    return $self->search_related_rs($annotation_source);
}

1;
