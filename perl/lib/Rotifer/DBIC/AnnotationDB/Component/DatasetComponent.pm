=head1 NAME

Rotifer::DBIC::AnnotationDB::Component::DatasetComponent - methods for datasets

=head1 DESCRIPTION

This modules extends Rotifer::DBIC::Result::Dataset.

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::DBIC::AnnotationDB::Component::DatasetComponent;

use strict;
use warnings;
use Carp::Clan qr/^Rotifer::DBIC::AnnotationDB/;
use Moose;
extends qw(
 DBIx::Class
);

=head2 add_related

 Title   : add_related
 Usage   : $text = $rs->add_related($something)
 Function: create relationships to row objects
 Returns : boolean
 Args    : row object (isa DBIx::Class::Row)

=cut

sub add_related {
    my ($self, $obj) = @_;
    return 0 unless (blessed $obj && $obj->isa("DBIx::Class::Row"));

    my $rel = undef;
    if (ref($obj) eq ref($self)) {
	$rel = $self->create_related("dataset_relationship_objects", { subject => $obj });
    } else {
	my $source = $obj->result_source;
	my $name   = $source->name;
	return undef unless $source->has_relationship("${name}_datasets");
	$rel = $self->create_related("${name}_datasets", { $name => $obj });
    }

    return $rel;
}

=head2 contains

 Title   : contains
 Usage   : $text = $rs->contains($something)
 Function: test if relationships to a row object exist
 Returns : boolean
 Args    : row object (isa DBIx::Class::Row)

=cut

sub contains {
    my ($self, $obj) = @_;
    return 0 unless (blessed $obj && $obj->isa("DBIx::Class::Row"));

    my $rel = undef;
    if (ref($obj) eq ref($self)) {
	$rel = $self->search_related_rs("dataset_relationship_objects", { subject => $obj });
    } else {
	my $source = $obj->result_source;
	my $name   = $source->name;
	return 0 unless $source->has_relationship("${name}_datasets");
	$rel = $self->search_related_rs("${name}_datasets", { $name => $obj });
    }

    return defined $rel && $rel->count;
}

__PACKAGE__->meta->make_immutable(inline_constructor => 0);
1;
