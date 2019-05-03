use utf8;
package Rotifer::DBIC::AnnotationDB::Result::Organism;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Rotifer::DBIC::AnnotationDB::Result::Organism

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

=item * L<Rotifer::DBIC::AnnotationDB::Component::Annotatable>

=item * L<Rotifer::DBIC::AnnotationDB::Component::WithDataset>

=back

=cut

__PACKAGE__->load_components(
  "Helper::ResultSet::SetOperations",
  "+Rotifer::DBIC::AnnotationDB::Component::Annotatable",
  "+Rotifer::DBIC::AnnotationDB::Component::WithDataset",
);

=head1 TABLE: C<organism>

=cut

__PACKAGE__->table("organism");

=head1 ACCESSORS

=head2 organism_id

  data_type: 'integer'
  is_auto_increment: 1
  is_nullable: 0

=head2 ncbi_taxon_id

  data_type: 'integer'
  is_nullable: 1

=head2 abbreviation

  data_type: 'varchar'
  is_nullable: 1
  size: 40

=head2 name

  data_type: 'varchar'
  is_nullable: 0
  size: 255

=head2 lineage

  data_type: 'text'
  is_nullable: 1

=head2 classification

  data_type: 'text'
  is_nullable: 1

=head2 comment

  data_type: 'text'
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "organism_id",
  { data_type => "integer", is_auto_increment => 1, is_nullable => 0 },
  "ncbi_taxon_id",
  { data_type => "integer", is_nullable => 1 },
  "abbreviation",
  { data_type => "varchar", is_nullable => 1, size => 40 },
  "name",
  { data_type => "varchar", is_nullable => 0, size => 255 },
  "lineage",
  { data_type => "text", is_nullable => 1 },
  "classification",
  { data_type => "text", is_nullable => 1 },
  "comment",
  { data_type => "text", is_nullable => 1 },
);

=head1 PRIMARY KEY

=over 4

=item * L</organism_id>

=back

=cut

__PACKAGE__->set_primary_key("organism_id");

=head1 UNIQUE CONSTRAINTS

=head2 C<organism_abbreviation_uniq>

=over 4

=item * L</abbreviation>

=back

=cut

__PACKAGE__->add_unique_constraint("organism_abbreviation_uniq", ["abbreviation"]);

=head2 C<organism_ncbitaxonid_uniq>

=over 4

=item * L</ncbi_taxon_id>

=back

=cut

__PACKAGE__->add_unique_constraint("organism_ncbitaxonid_uniq", ["ncbi_taxon_id"]);

=head1 RELATIONS

=head2 biodatas

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::Biodata>

=cut

__PACKAGE__->has_many(
  "biodatas",
  "Rotifer::DBIC::AnnotationDB::Result::Biodata",
  { "foreign.organism_id" => "self.organism_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 organism_attribute_values

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::OrganismAttributeValue>

=cut

__PACKAGE__->has_many(
  "organism_attribute_values",
  "Rotifer::DBIC::AnnotationDB::Result::OrganismAttributeValue",
  { "foreign.organism_id" => "self.organism_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 organism_datasets

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::OrganismDataset>

=cut

__PACKAGE__->has_many(
  "organism_datasets",
  "Rotifer::DBIC::AnnotationDB::Result::OrganismDataset",
  { "foreign.organism_id" => "self.organism_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 organism_dbxrefs

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::OrganismDbxref>

=cut

__PACKAGE__->has_many(
  "organism_dbxrefs",
  "Rotifer::DBIC::AnnotationDB::Result::OrganismDbxref",
  { "foreign.organism_id" => "self.organism_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);


# Created by DBIx::Class::Schema::Loader v0.07033 @ 2013-06-17 13:33:01
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:AsAKInaF05Jafz8zWzI09Q


# You can replace this text with custom code or comments, and it will be preserved on regeneration
__PACKAGE__->meta->make_immutable;
1;
