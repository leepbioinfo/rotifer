use utf8;
package Rotifer::DBIC::AnnotationDB::Result::Biodata;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Rotifer::DBIC::AnnotationDB::Result::Biodata

=cut

use strict;
use warnings;

use Moose;
use MooseX::NonMoose;
use MooseX::MarkAsMethods autoclean => 1;
extends 'DBIx::Class::Core';

=head1 COMPONENTS LOADED

=over 4

=item * L<Rotifer::DBIC::AnnotationDB::Component::DynamicSubclassByTerm>

=item * L<Rotifer::DBIC::AnnotationDB::Component::BiodataComponent>

=item * L<Rotifer::DBIC::AnnotationDB::Component::Annotatable>

=item * L<DBIx::Class::Helper::ResultSet::SetOperations>

=item * L<DBIx::Class::InflateColumn::DateTime>

=item * L<DBIx::Class::Core>

=back

=cut

__PACKAGE__->load_components(
  "+Rotifer::DBIC::AnnotationDB::Component::DynamicSubclassByTerm",
  "+Rotifer::DBIC::AnnotationDB::Component::BiodataComponent",
  "+Rotifer::DBIC::AnnotationDB::Component::Annotatable",
  "Helper::ResultSet::SetOperations",
  "InflateColumn::DateTime",
  "Core",
);

=head1 TABLE: C<biodata>

=cut

__PACKAGE__->table("biodata");

=head1 ACCESSORS

=head2 biodata_id

  data_type: 'integer'
  is_auto_increment: 1
  is_nullable: 0

=head2 term_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=head2 organism_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 1

=head2 biosequence_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 1

=head2 identifier

  data_type: 'varchar'
  is_nullable: 0
  size: 255

=head2 rank

  data_type: 'smallint'
  default_value: 0
  is_nullable: 0

=head2 is_rawsequence

  data_type: 'tinyint'
  default_value: 0
  is_nullable: 0

=head2 is_obsolete

  data_type: 'tinyint'
  default_value: 0
  is_nullable: 0

=head2 description

  data_type: 'text'
  is_nullable: 1

=head2 lastmodified

  data_type: 'timestamp'
  datetime_undef_if_invalid: 1
  default_value: current_timestamp
  is_nullable: 0

=cut

__PACKAGE__->add_columns(
  "biodata_id",
  { data_type => "integer", is_auto_increment => 1, is_nullable => 0 },
  "term_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
  "organism_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 1 },
  "biosequence_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 1 },
  "identifier",
  { data_type => "varchar", is_nullable => 0, size => 255 },
  "rank",
  { data_type => "smallint", default_value => 0, is_nullable => 0 },
  "is_rawsequence",
  { data_type => "tinyint", default_value => 0, is_nullable => 0 },
  "is_obsolete",
  { data_type => "tinyint", default_value => 0, is_nullable => 0 },
  "description",
  { data_type => "text", is_nullable => 1 },
  "lastmodified",
  {
    data_type => "timestamp",
    datetime_undef_if_invalid => 1,
    default_value => \"current_timestamp",
    is_nullable => 0,
  },
);

=head1 PRIMARY KEY

=over 4

=item * L</biodata_id>

=back

=cut

__PACKAGE__->set_primary_key("biodata_id");

=head1 UNIQUE CONSTRAINTS

=head2 C<biodata_uniq>

=over 4

=item * L</term_id>

=item * L</identifier>

=item * L</rank>

=back

=cut

__PACKAGE__->add_unique_constraint("biodata_uniq", ["term_id", "identifier", "rank"]);

=head1 RELATIONS

=head2 biodata_attribute_values

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::BiodataAttributeValue>

=cut

__PACKAGE__->has_many(
  "biodata_attribute_values",
  "Rotifer::DBIC::AnnotationDB::Result::BiodataAttributeValue",
  { "foreign.biodata_id" => "self.biodata_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 biodata_datasets

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::BiodataDataset>

=cut

__PACKAGE__->has_many(
  "biodata_datasets",
  "Rotifer::DBIC::AnnotationDB::Result::BiodataDataset",
  { "foreign.biodata_id" => "self.biodata_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 biodata_dbxrefs

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::BiodataDbxref>

=cut

__PACKAGE__->has_many(
  "biodata_dbxrefs",
  "Rotifer::DBIC::AnnotationDB::Result::BiodataDbxref",
  { "foreign.biodata_id" => "self.biodata_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 biodata_relationship_objects

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::BiodataRelationship>

=cut

__PACKAGE__->has_many(
  "biodata_relationship_objects",
  "Rotifer::DBIC::AnnotationDB::Result::BiodataRelationship",
  { "foreign.object_id" => "self.biodata_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 biodata_relationship_subjects

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::BiodataRelationship>

=cut

__PACKAGE__->has_many(
  "biodata_relationship_subjects",
  "Rotifer::DBIC::AnnotationDB::Result::BiodataRelationship",
  { "foreign.subject_id" => "self.biodata_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 biosequence

Type: belongs_to

Related object: L<Rotifer::DBIC::AnnotationDB::Result::Biosequence>

=cut

__PACKAGE__->belongs_to(
  "biosequence",
  "Rotifer::DBIC::AnnotationDB::Result::Biosequence",
  { biosequence_id => "biosequence_id" },
  {
    is_deferrable => 1,
    join_type     => "LEFT",
    on_delete     => "CASCADE",
    on_update     => "CASCADE",
  },
);

=head2 organism

Type: belongs_to

Related object: L<Rotifer::DBIC::AnnotationDB::Result::Organism>

=cut

__PACKAGE__->belongs_to(
  "organism",
  "Rotifer::DBIC::AnnotationDB::Result::Organism",
  { organism_id => "organism_id" },
  {
    is_deferrable => 1,
    join_type     => "LEFT",
    on_delete     => "CASCADE",
    on_update     => "CASCADE",
  },
);

=head2 term

Type: belongs_to

Related object: L<Rotifer::DBIC::AnnotationDB::Result::Term>

=cut

__PACKAGE__->belongs_to(
  "term",
  "Rotifer::DBIC::AnnotationDB::Result::Term",
  { term_id => "term_id" },
  { is_deferrable => 1, on_delete => "NO ACTION", on_update => "NO ACTION" },
);

=head2 datasets

Type: many_to_many

Composing rels: L</biodata_datasets> -> dataset

=cut

__PACKAGE__->many_to_many("datasets", "biodata_datasets", "dataset");

=head2 dbxrefs

Type: many_to_many

Composing rels: L</biodata_dbxrefs> -> dbxref

=cut

__PACKAGE__->many_to_many("dbxrefs", "biodata_dbxrefs", "dbxref");


# Created by DBIx::Class::Schema::Loader v0.07033 @ 2013-05-22 17:01:22
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:N7KX/3uCUGNn9Qa6ygrHlw


__PACKAGE__->meta->make_immutable;
1;
