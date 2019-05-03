use utf8;
package Rotifer::DBIC::AnnotationDB::Result::Dbxref;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Rotifer::DBIC::AnnotationDB::Result::Dbxref

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

=item * L<Rotifer::DBIC::AnnotationDB::Component::WithDataset>

=back

=cut

__PACKAGE__->load_components(
  "Helper::ResultSet::SetOperations",
  "+Rotifer::DBIC::AnnotationDB::Component::WithDataset",
);

=head1 TABLE: C<dbxref>

=cut

__PACKAGE__->table("dbxref");

=head1 ACCESSORS

=head2 dbxref_id

  data_type: 'integer'
  is_auto_increment: 1
  is_nullable: 0

=head2 dbname

  data_type: 'varchar'
  is_nullable: 0
  size: 40

=head2 accession

  data_type: 'varchar'
  is_nullable: 0
  size: 128

=head2 version

  data_type: 'smallint'
  default_value: 0
  is_nullable: 0

=head2 acctype

  data_type: 'varchar'
  is_nullable: 1
  size: 40

=head2 accgroup

  data_type: 'varchar'
  is_nullable: 1
  size: 128

=cut

__PACKAGE__->add_columns(
  "dbxref_id",
  { data_type => "integer", is_auto_increment => 1, is_nullable => 0 },
  "dbname",
  { data_type => "varchar", is_nullable => 0, size => 40 },
  "accession",
  { data_type => "varchar", is_nullable => 0, size => 128 },
  "version",
  { data_type => "smallint", default_value => 0, is_nullable => 0 },
  "acctype",
  { data_type => "varchar", is_nullable => 1, size => 40 },
  "accgroup",
  { data_type => "varchar", is_nullable => 1, size => 128 },
);

=head1 PRIMARY KEY

=over 4

=item * L</dbxref_id>

=back

=cut

__PACKAGE__->set_primary_key("dbxref_id");

=head1 UNIQUE CONSTRAINTS

=head2 C<dbxref_uniq>

=over 4

=item * L</accession>

=item * L</dbname>

=item * L</version>

=back

=cut

__PACKAGE__->add_unique_constraint("dbxref_uniq", ["accession", "dbname", "version"]);

=head1 RELATIONS

=head2 biodata_dbxrefs

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::BiodataDbxref>

=cut

__PACKAGE__->has_many(
  "biodata_dbxrefs",
  "Rotifer::DBIC::AnnotationDB::Result::BiodataDbxref",
  { "foreign.dbxref_id" => "self.dbxref_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 biosequence_dbxrefs

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::BiosequenceDbxref>

=cut

__PACKAGE__->has_many(
  "biosequence_dbxrefs",
  "Rotifer::DBIC::AnnotationDB::Result::BiosequenceDbxref",
  { "foreign.dbxref_id" => "self.dbxref_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 dbxref_datasets

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::DbxrefDataset>

=cut

__PACKAGE__->has_many(
  "dbxref_datasets",
  "Rotifer::DBIC::AnnotationDB::Result::DbxrefDataset",
  { "foreign.dbxref_id" => "self.dbxref_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 dbxref_qualifier_values

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::DbxrefQualifierValue>

=cut

__PACKAGE__->has_many(
  "dbxref_qualifier_values",
  "Rotifer::DBIC::AnnotationDB::Result::DbxrefQualifierValue",
  { "foreign.dbxref_id" => "self.dbxref_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 organism_dbxrefs

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::OrganismDbxref>

=cut

__PACKAGE__->has_many(
  "organism_dbxrefs",
  "Rotifer::DBIC::AnnotationDB::Result::OrganismDbxref",
  { "foreign.dbxref_id" => "self.dbxref_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 term_dbxrefs

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::TermDbxref>

=cut

__PACKAGE__->has_many(
  "term_dbxrefs",
  "Rotifer::DBIC::AnnotationDB::Result::TermDbxref",
  { "foreign.dbxref_id" => "self.dbxref_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 biodatas

Type: many_to_many

Composing rels: L</biodata_dbxrefs> -> biodata

=cut

__PACKAGE__->many_to_many("biodatas", "biodata_dbxrefs", "biodata");

=head2 biosequences

Type: many_to_many

Composing rels: L</biosequence_dbxrefs> -> biosequence

=cut

__PACKAGE__->many_to_many("biosequences", "biosequence_dbxrefs", "biosequence");

=head2 datasets

Type: many_to_many

Composing rels: L</dbxref_datasets> -> dataset

=cut

__PACKAGE__->many_to_many("datasets", "dbxref_datasets", "dataset");


# Created by DBIx::Class::Schema::Loader v0.07033 @ 2013-06-17 13:33:01
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:LHx3Zj4mCGN2ypAypR4K0g


# You can replace this text with custom code or comments, and it will be preserved on regeneration
__PACKAGE__->meta->make_immutable;
1;
