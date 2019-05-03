use utf8;
package Rotifer::DBIC::NCBI_GENE::Result::GeneGroup;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Rotifer::DBIC::NCBI_GENE::Result::GeneGroup

=cut

use strict;
use warnings;

use Moose;
use MooseX::NonMoose;
use MooseX::MarkAsMethods autoclean => 1;
extends 'DBIx::Class::Core';

=head1 TABLE: C<gene_group>

=cut

__PACKAGE__->table("gene_group");

=head1 ACCESSORS

=head2 ncbi_taxon_id

  data_type: 'integer'
  is_nullable: 0

=head2 gene_id

  data_type: 'integer'
  is_nullable: 0

=head2 relationship

  data_type: 'text'
  is_nullable: 1

=head2 other_taxon_id

  data_type: 'integer'
  is_nullable: 0

=head2 other_gene_id

  data_type: 'integer'
  is_nullable: 0

=cut

__PACKAGE__->add_columns(
  "ncbi_taxon_id",
  { data_type => "integer", is_nullable => 0 },
  "gene_id",
  { data_type => "integer", is_nullable => 0 },
  "relationship",
  { data_type => "text", is_nullable => 1 },
  "other_taxon_id",
  { data_type => "integer", is_nullable => 0 },
  "other_gene_id",
  { data_type => "integer", is_nullable => 0 },
);


# Created by DBIx::Class::Schema::Loader v0.07035 @ 2013-07-18 14:44:41
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:Dg2TesZxgW5Ybj5xVzI8+Q


# You can replace this text with custom code or comments, and it will be preserved on regeneration
__PACKAGE__->meta->make_immutable;
1;
