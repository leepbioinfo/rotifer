use utf8;
package Rotifer::DBIC::AnnotationDB::Result::BiodataRelationship;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Rotifer::DBIC::AnnotationDB::Result::BiodataRelationship

=cut

use strict;
use warnings;

use Moose;
use MooseX::NonMoose;
use MooseX::MarkAsMethods autoclean => 1;
extends 'DBIx::Class::Core';

=head1 COMPONENTS LOADED

=over 4

=item * L<Rotifer::DBIC::AnnotationDB::Component::BiodataRelationshipComponent>

=item * L<Rotifer::DBIC::AnnotationDB::Component::Annotatable>

=item * L<DBIx::Class::Helper::ResultSet::SetOperations>

=back

=cut

__PACKAGE__->load_components(
  "+Rotifer::DBIC::AnnotationDB::Component::BiodataRelationshipComponent",
  "+Rotifer::DBIC::AnnotationDB::Component::Annotatable",
  "Helper::ResultSet::SetOperations",
);

=head1 TABLE: C<biodata_relationship>

=cut

__PACKAGE__->table("biodata_relationship");

=head1 ACCESSORS

=head2 biodata_relationship_id

  data_type: 'integer'
  is_auto_increment: 1
  is_nullable: 0

=head2 subject_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=head2 predicate_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=head2 object_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=head2 rank

  data_type: 'integer'
  default_value: 0
  is_nullable: 0

=head2 start_coord

  data_type: 'integer'
  is_nullable: 1

=head2 end_coord

  data_type: 'integer'
  is_nullable: 1

=head2 strand

  data_type: 'smallint'
  default_value: 0
  is_nullable: 1

=head2 score

  data_type: 'double precision'
  default_value: 0
  is_nullable: 1

=head2 is_location

  data_type: 'tinyint'
  default_value: 0
  is_nullable: 0

=cut

__PACKAGE__->add_columns(
  "biodata_relationship_id",
  { data_type => "integer", is_auto_increment => 1, is_nullable => 0 },
  "subject_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
  "predicate_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
  "object_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
  "rank",
  { data_type => "integer", default_value => 0, is_nullable => 0 },
  "start_coord",
  { data_type => "integer", is_nullable => 1 },
  "end_coord",
  { data_type => "integer", is_nullable => 1 },
  "strand",
  { data_type => "smallint", default_value => 0, is_nullable => 1 },
  "score",
  { data_type => "double precision", default_value => 0, is_nullable => 1 },
  "is_location",
  { data_type => "tinyint", default_value => 0, is_nullable => 0 },
);

=head1 PRIMARY KEY

=over 4

=item * L</biodata_relationship_id>

=back

=cut

__PACKAGE__->set_primary_key("biodata_relationship_id");

=head1 UNIQUE CONSTRAINTS

=head2 C<biodata_relationship_uniq>

=over 4

=item * L</object_id>

=item * L</subject_id>

=item * L</predicate_id>

=item * L</rank>

=back

=cut

__PACKAGE__->add_unique_constraint(
  "biodata_relationship_uniq",
  ["object_id", "subject_id", "predicate_id", "rank"],
);

=head1 RELATIONS

=head2 biodata_relationship_attribute_values

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::BiodataRelationshipAttributeValue>

=cut

__PACKAGE__->has_many(
  "biodata_relationship_attribute_values",
  "Rotifer::DBIC::AnnotationDB::Result::BiodataRelationshipAttributeValue",
  {
    "foreign.biodata_relationship_id" => "self.biodata_relationship_id",
  },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 biodata_relationship_datasets

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::BiodataRelationshipDataset>

=cut

__PACKAGE__->has_many(
  "biodata_relationship_datasets",
  "Rotifer::DBIC::AnnotationDB::Result::BiodataRelationshipDataset",
  {
    "foreign.biodata_relationship_id" => "self.biodata_relationship_id",
  },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 object

Type: belongs_to

Related object: L<Rotifer::DBIC::AnnotationDB::Result::Biodata>

=cut

__PACKAGE__->belongs_to(
  "object",
  "Rotifer::DBIC::AnnotationDB::Result::Biodata",
  { biodata_id => "object_id" },
  { is_deferrable => 1, on_delete => "CASCADE", on_update => "CASCADE" },
);

=head2 predicate

Type: belongs_to

Related object: L<Rotifer::DBIC::AnnotationDB::Result::Term>

=cut

__PACKAGE__->belongs_to(
  "predicate",
  "Rotifer::DBIC::AnnotationDB::Result::Term",
  { term_id => "predicate_id" },
  { is_deferrable => 1, on_delete => "CASCADE", on_update => "CASCADE" },
);

=head2 subject

Type: belongs_to

Related object: L<Rotifer::DBIC::AnnotationDB::Result::Biodata>

=cut

__PACKAGE__->belongs_to(
  "subject",
  "Rotifer::DBIC::AnnotationDB::Result::Biodata",
  { biodata_id => "subject_id" },
  { is_deferrable => 1, on_delete => "CASCADE", on_update => "CASCADE" },
);

=head2 datasets

Type: many_to_many

Composing rels: L</biodata_relationship_datasets> -> dataset

=cut

__PACKAGE__->many_to_many("datasets", "biodata_relationship_datasets", "dataset");


# Created by DBIx::Class::Schema::Loader v0.07033 @ 2013-06-17 01:39:13
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:gYVYGihl9T56NxfwY/lXjA


# You can replace this text with custom code or comments, and it will be preserved on regeneration
__PACKAGE__->meta->make_immutable;
1;
