use utf8;
package Rotifer::DBIC::AnnotationDB::Result::Dataset;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Rotifer::DBIC::AnnotationDB::Result::Dataset

=cut

use strict;
use warnings;

use Moose;
use MooseX::NonMoose;
use MooseX::MarkAsMethods autoclean => 1;
extends 'DBIx::Class::Core';

=head1 COMPONENTS LOADED

=over 4

=item * L<DBIx::Class::Helper::ResultSet::SetOperations>

=item * L<DBIx::Class::InflateColumn::DateTime>

=item * L<Rotifer::DBIC::AnnotationDB::Component::DatasetComponent>

=item * L<Rotifer::DBIC::AnnotationDB::Component::Annotatable>

=back

=cut

__PACKAGE__->load_components(
  "Helper::ResultSet::SetOperations",
  "InflateColumn::DateTime",
  "+Rotifer::DBIC::AnnotationDB::Component::DatasetComponent",
  "+Rotifer::DBIC::AnnotationDB::Component::Annotatable",
);

=head1 TABLE: C<dataset>

=cut

__PACKAGE__->table("dataset");

=head1 ACCESSORS

=head2 dataset_id

  data_type: 'integer'
  is_auto_increment: 1
  is_nullable: 0

=head2 name

  data_type: 'varchar'
  is_nullable: 0
  size: 255

=head2 authority

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 lastmodified

  data_type: 'timestamp'
  datetime_undef_if_invalid: 1
  default_value: current_timestamp
  is_nullable: 0

=head2 loaded_by

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 loaded_from

  data_type: 'text'
  is_nullable: 1

=head2 description

  data_type: 'text'
  is_nullable: 1

=head2 is_project

  data_type: 'tinyint'
  default_value: 0
  is_nullable: 0

=cut

__PACKAGE__->add_columns(
  "dataset_id",
  { data_type => "integer", is_auto_increment => 1, is_nullable => 0 },
  "name",
  { data_type => "varchar", is_nullable => 0, size => 255 },
  "authority",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "lastmodified",
  {
    data_type => "timestamp",
    datetime_undef_if_invalid => 1,
    default_value => \"current_timestamp",
    is_nullable => 0,
  },
  "loaded_by",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "loaded_from",
  { data_type => "text", is_nullable => 1 },
  "description",
  { data_type => "text", is_nullable => 1 },
  "is_project",
  { data_type => "tinyint", default_value => 0, is_nullable => 0 },
);

=head1 PRIMARY KEY

=over 4

=item * L</dataset_id>

=back

=cut

__PACKAGE__->set_primary_key("dataset_id");

=head1 UNIQUE CONSTRAINTS

=head2 C<dataset_uniq>

=over 4

=item * L</name>

=back

=cut

__PACKAGE__->add_unique_constraint("dataset_uniq", ["name"]);

=head1 RELATIONS

=head2 biodata_datasets

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::BiodataDataset>

=cut

__PACKAGE__->has_many(
  "biodata_datasets",
  "Rotifer::DBIC::AnnotationDB::Result::BiodataDataset",
  { "foreign.dataset_id" => "self.dataset_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 biodata_relationship_datasets

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::BiodataRelationshipDataset>

=cut

__PACKAGE__->has_many(
  "biodata_relationship_datasets",
  "Rotifer::DBIC::AnnotationDB::Result::BiodataRelationshipDataset",
  { "foreign.dataset_id" => "self.dataset_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 biosequence_datasets

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::BiosequenceDataset>

=cut

__PACKAGE__->has_many(
  "biosequence_datasets",
  "Rotifer::DBIC::AnnotationDB::Result::BiosequenceDataset",
  { "foreign.dataset_id" => "self.dataset_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 dataset_attribute_values

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::DatasetAttributeValue>

=cut

__PACKAGE__->has_many(
  "dataset_attribute_values",
  "Rotifer::DBIC::AnnotationDB::Result::DatasetAttributeValue",
  { "foreign.dataset_id" => "self.dataset_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 dataset_relationship_objects

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::DatasetRelationship>

=cut

__PACKAGE__->has_many(
  "dataset_relationship_objects",
  "Rotifer::DBIC::AnnotationDB::Result::DatasetRelationship",
  { "foreign.object_id" => "self.dataset_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 dataset_relationship_subjects

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::DatasetRelationship>

=cut

__PACKAGE__->has_many(
  "dataset_relationship_subjects",
  "Rotifer::DBIC::AnnotationDB::Result::DatasetRelationship",
  { "foreign.subject_id" => "self.dataset_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 dbxref_datasets

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::DbxrefDataset>

=cut

__PACKAGE__->has_many(
  "dbxref_datasets",
  "Rotifer::DBIC::AnnotationDB::Result::DbxrefDataset",
  { "foreign.dataset_id" => "self.dataset_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 organism_datasets

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::OrganismDataset>

=cut

__PACKAGE__->has_many(
  "organism_datasets",
  "Rotifer::DBIC::AnnotationDB::Result::OrganismDataset",
  { "foreign.dataset_id" => "self.dataset_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 biodata_relationships

Type: many_to_many

Composing rels: L</biodata_relationship_datasets> -> biodata_relationship

=cut

__PACKAGE__->many_to_many(
  "biodata_relationships",
  "biodata_relationship_datasets",
  "biodata_relationship",
);

=head2 biodatas

Type: many_to_many

Composing rels: L</biodata_datasets> -> biodata

=cut

__PACKAGE__->many_to_many("biodatas", "biodata_datasets", "biodata");

=head2 biosequences

Type: many_to_many

Composing rels: L</biosequence_datasets> -> biosequence

=cut

__PACKAGE__->many_to_many("biosequences", "biosequence_datasets", "biosequence");

=head2 dbxrefs

Type: many_to_many

Composing rels: L</dbxref_datasets> -> dbxref

=cut

__PACKAGE__->many_to_many("dbxrefs", "dbxref_datasets", "dbxref");


# Created by DBIx::Class::Schema::Loader v0.07033 @ 2013-06-17 13:33:01
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:Xm/K5fDBUOfRQOIfKmDCIg


# You can replace this text with custom code or comments, and it will be preserved on regeneration
__PACKAGE__->meta->make_immutable;
1;
