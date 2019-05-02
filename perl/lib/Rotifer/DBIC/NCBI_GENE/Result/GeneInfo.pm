use utf8;
package Rotifer::DBIC::NCBI_GENE::Result::GeneInfo;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Rotifer::DBIC::NCBI_GENE::Result::GeneInfo

=cut

use strict;
use warnings;

use Moose;
use MooseX::NonMoose;
use MooseX::MarkAsMethods autoclean => 1;
extends 'DBIx::Class::Core';

=head1 TABLE: C<gene_info>

=cut

__PACKAGE__->table("gene_info");

=head1 ACCESSORS

=head2 ncbi_taxon_id

  data_type: 'integer'
  is_nullable: 0

=head2 gene_id

  data_type: 'integer'
  is_nullable: 0

=head2 symbol

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 locus_tag

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 synonyms

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 dbxrefs

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 chromosome

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 map_location

  data_type: 'text'
  is_nullable: 1

=head2 description

  data_type: 'text'
  is_nullable: 1

=head2 type_of_gene

  data_type: 'varchar'
  is_nullable: 1
  size: 64

=head2 symbol_from_authority

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 full_name_from_authority

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 nomenclature_status

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 other_designations

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 modification_date

  data_type: 'integer'
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "ncbi_taxon_id",
  { data_type => "integer", is_nullable => 0 },
  "gene_id",
  { data_type => "integer", is_nullable => 0 },
  "symbol",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "locus_tag",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "synonyms",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "dbxrefs",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "chromosome",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "map_location",
  { data_type => "text", is_nullable => 1 },
  "description",
  { data_type => "text", is_nullable => 1 },
  "type_of_gene",
  { data_type => "varchar", is_nullable => 1, size => 64 },
  "symbol_from_authority",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "full_name_from_authority",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "nomenclature_status",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "other_designations",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "modification_date",
  { data_type => "integer", is_nullable => 1 },
);

=head1 PRIMARY KEY

=over 4

=item * L</gene_id>

=back

=cut

__PACKAGE__->set_primary_key("gene_id");


# Created by DBIx::Class::Schema::Loader v0.07035 @ 2013-07-18 14:44:41
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:D6APa2A2IPb+SBc1nlDc8w


# You can replace this text with custom code or comments, and it will be preserved on regeneration
__PACKAGE__->meta->make_immutable;
1;
