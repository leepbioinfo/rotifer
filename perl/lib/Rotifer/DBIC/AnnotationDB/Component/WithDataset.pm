=head1 NAME

Rotifer::DBIC::AnnotationDB::Component::WithDataset - 
 automatically handle relationships with datasets

=head1 DESCRIPTION

This modules takes care of the relationships between rows in a table
and the current dataset.

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::DBIC::AnnotationDB::Component::WithDataset;

use strict;
use warnings;
use Carp::Clan qr/^Rotifer::DBIC::AnnotationDB/;
use Data::UUID;
use base qw(DBIx::Class);

sub insert {
    my ( $class, $attrs ) = @_;

    # Create object
    my $obj = $class->next::method($attrs);

    # Automatically add to datasets
    my $relname = $obj->result_source->name . "_datasets";
    if ($obj->result_source->has_relationship($relname)) {
	my $schema = $obj->result_source->schema;
	foreach my $set (@{ $schema->datasets }) {
	    $obj->find_or_create_related($relname, { dataset => $set })
	};
    }

    return $obj;
}

1;
