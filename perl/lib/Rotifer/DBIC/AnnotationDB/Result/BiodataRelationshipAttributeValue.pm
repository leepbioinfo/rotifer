use utf8;
package Rotifer::DBIC::AnnotationDB::Result::BiodataRelationshipAttributeValue;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Rotifer::DBIC::AnnotationDB::Result::BiodataRelationshipAttributeValue

=cut

use strict;
use warnings;

use Moose;
use MooseX::NonMoose;
use MooseX::MarkAsMethods autoclean => 1;
extends 'DBIx::Class::Core';

=head1 TABLE: C<biodata_relationship_attribute_value>

=cut

__PACKAGE__->table("biodata_relationship_attribute_value");

=head1 ACCESSORS

=head2 biodata_relationship_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=head2 term_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=head2 rank

  data_type: 'integer'
  default_value: 0
  is_nullable: 0

=head2 value

  data_type: 'text'
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "biodata_relationship_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
  "term_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
  "rank",
  { data_type => "integer", default_value => 0, is_nullable => 0 },
  "value",
  { data_type => "text", is_nullable => 1 },
);

=head1 PRIMARY KEY

=over 4

=item * L</biodata_relationship_id>

=item * L</term_id>

=item * L</rank>

=back

=cut

__PACKAGE__->set_primary_key("biodata_relationship_id", "term_id", "rank");

=head1 RELATIONS

=head2 biodata_relationship

Type: belongs_to

Related object: L<Rotifer::DBIC::AnnotationDB::Result::BiodataRelationship>

=cut

__PACKAGE__->belongs_to(
  "biodata_relationship",
  "Rotifer::DBIC::AnnotationDB::Result::BiodataRelationship",
  { biodata_relationship_id => "biodata_relationship_id" },
  { is_deferrable => 1, on_delete => "NO ACTION", on_update => "NO ACTION" },
);

=head2 term

Type: belongs_to

Related object: L<Rotifer::DBIC::AnnotationDB::Result::Term>

=cut

__PACKAGE__->belongs_to(
  "term",
  "Rotifer::DBIC::AnnotationDB::Result::Term",
  { term_id => "term_id" },
  { is_deferrable => 1, on_delete => "CASCADE", on_update => "CASCADE" },
);


# Created by DBIx::Class::Schema::Loader v0.07033 @ 2013-05-08 16:20:32
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:X1YVoNMiZW5YIkqgU1l77Q


# You can replace this text with custom code or comments, and it will be preserved on regeneration
__PACKAGE__->meta->make_immutable;
1;
