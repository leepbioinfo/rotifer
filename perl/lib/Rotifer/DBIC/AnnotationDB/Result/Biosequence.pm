use utf8;
package Rotifer::DBIC::AnnotationDB::Result::Biosequence;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Rotifer::DBIC::AnnotationDB::Result::Biosequence

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

=item * L<Rotifer::DBIC::AnnotationDB::Component::BiosequenceComponent>

=back

=cut

__PACKAGE__->load_components(
  "Helper::ResultSet::SetOperations",
  "+Rotifer::DBIC::AnnotationDB::Component::BiosequenceComponent",
);

=head1 TABLE: C<biosequence>

=cut

__PACKAGE__->table("biosequence");

=head1 ACCESSORS

=head2 biosequence_id

  data_type: 'integer'
  is_auto_increment: 1
  is_nullable: 0

=head2 term_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=head2 checksum

  data_type: 'varchar'
  is_nullable: 0
  size: 255

=head2 seqlen

  data_type: 'integer'
  is_nullable: 0

=head2 circular

  data_type: 'tinyint'
  default_value: 0
  is_nullable: 0

=head2 residues

  data_type: 'text'
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "biosequence_id",
  { data_type => "integer", is_auto_increment => 1, is_nullable => 0 },
  "term_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
  "checksum",
  { data_type => "varchar", is_nullable => 0, size => 255 },
  "seqlen",
  { data_type => "integer", is_nullable => 0 },
  "circular",
  { data_type => "tinyint", default_value => 0, is_nullable => 0 },
  "residues",
  { data_type => "text", is_nullable => 1 },
);

=head1 PRIMARY KEY

=over 4

=item * L</biosequence_id>

=back

=cut

__PACKAGE__->set_primary_key("biosequence_id");

=head1 UNIQUE CONSTRAINTS

=head2 C<biosequence_uniq>

=over 4

=item * L</term_id>

=item * L</checksum>

=back

=cut

__PACKAGE__->add_unique_constraint("biosequence_uniq", ["term_id", "checksum"]);

=head1 RELATIONS

=head2 biodatas

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::Biodata>

=cut

__PACKAGE__->has_many(
  "biodatas",
  "Rotifer::DBIC::AnnotationDB::Result::Biodata",
  { "foreign.biosequence_id" => "self.biosequence_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 biosequence_datasets

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::BiosequenceDataset>

=cut

__PACKAGE__->has_many(
  "biosequence_datasets",
  "Rotifer::DBIC::AnnotationDB::Result::BiosequenceDataset",
  { "foreign.biosequence_id" => "self.biosequence_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 biosequence_dbxrefs

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::BiosequenceDbxref>

=cut

__PACKAGE__->has_many(
  "biosequence_dbxrefs",
  "Rotifer::DBIC::AnnotationDB::Result::BiosequenceDbxref",
  { "foreign.biosequence_id" => "self.biosequence_id" },
  { cascade_copy => 0, cascade_delete => 0 },
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

Composing rels: L</biosequence_datasets> -> dataset

=cut

__PACKAGE__->many_to_many("datasets", "biosequence_datasets", "dataset");

=head2 dbxrefs

Type: many_to_many

Composing rels: L</biosequence_dbxrefs> -> dbxref

=cut

__PACKAGE__->many_to_many("dbxrefs", "biosequence_dbxrefs", "dbxref");


# Created by DBIx::Class::Schema::Loader v0.07033 @ 2013-06-17 01:39:13
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:Ja8HKspBfUhxRZnCRualFA


# You can replace this text with custom code or comments, and it will be preserved on regeneration
__PACKAGE__->meta->make_immutable;
1;
