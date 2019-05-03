use utf8;
package Rotifer::DBIC::NCBI_GENE::Result::GeneNeighbor;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Rotifer::DBIC::NCBI_GENE::Result::GeneNeighbor

=cut

use strict;
use warnings;

use Moose;
use MooseX::NonMoose;
use MooseX::MarkAsMethods autoclean => 1;
extends 'DBIx::Class::Core';

=head1 TABLE: C<gene_neighbors>

=cut

__PACKAGE__->table("gene_neighbors");

=head1 ACCESSORS

=head2 ncbi_taxon_id

  data_type: 'integer'
  is_nullable: 0

=head2 gene_id

  data_type: 'integer'
  is_nullable: 0

=head2 genomic_accession

  data_type: 'varchar'
  is_nullable: 1
  size: 64

=head2 genomic_gi

  data_type: 'integer'
  is_nullable: 1

=head2 start_pos

  data_type: 'integer'
  is_nullable: 1

=head2 end_pos

  data_type: 'integer'
  is_nullable: 1

=head2 orientation

  data_type: 'char'
  is_nullable: 1
  size: 1

=head2 chromosome

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 left_geneids

  data_type: 'text'
  is_nullable: 1

=head2 distance_to_left

  data_type: 'integer'
  is_nullable: 1

=head2 right_geneids

  data_type: 'text'
  is_nullable: 1

=head2 distance_to_right

  data_type: 'integer'
  is_nullable: 1

=head2 overlapping_geneids

  data_type: 'text'
  is_nullable: 1

=head2 assembly

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=cut

__PACKAGE__->add_columns(
  "ncbi_taxon_id",
  { data_type => "integer", is_nullable => 0 },
  "gene_id",
  { data_type => "integer", is_nullable => 0 },
  "genomic_accession",
  { data_type => "varchar", is_nullable => 1, size => 64 },
  "genomic_gi",
  { data_type => "integer", is_nullable => 1 },
  "start_pos",
  { data_type => "integer", is_nullable => 1 },
  "end_pos",
  { data_type => "integer", is_nullable => 1 },
  "orientation",
  { data_type => "char", is_nullable => 1, size => 1 },
  "chromosome",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "left_geneids",
  { data_type => "text", is_nullable => 1 },
  "distance_to_left",
  { data_type => "integer", is_nullable => 1 },
  "right_geneids",
  { data_type => "text", is_nullable => 1 },
  "distance_to_right",
  { data_type => "integer", is_nullable => 1 },
  "overlapping_geneids",
  { data_type => "text", is_nullable => 1 },
  "assembly",
  { data_type => "varchar", is_nullable => 1, size => 255 },
);


# Created by DBIx::Class::Schema::Loader v0.07035 @ 2013-07-18 14:44:41
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:oW+zryG0eQU/UFVAowEkwA


# You can replace this text with custom code or comments, and it will be preserved on regeneration
__PACKAGE__->meta->make_immutable;
1;
