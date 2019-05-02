use utf8;
package Rotifer::DBIC::NCBI_GENE::Result::GeneHistory;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Rotifer::DBIC::NCBI_GENE::Result::GeneHistory

=cut

use strict;
use warnings;

use Moose;
use MooseX::NonMoose;
use MooseX::MarkAsMethods autoclean => 1;
extends 'DBIx::Class::Core';

=head1 TABLE: C<gene_history>

=cut

__PACKAGE__->table("gene_history");

=head1 ACCESSORS

=head2 ncbi_taxon_id

  data_type: 'integer'
  is_nullable: 0

=head2 gene_id

  data_type: 'integer'
  is_nullable: 0

=head2 discontinued_gene_id

  data_type: 'integer'
  is_nullable: 1

=head2 discontinued_symbol

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 discontinue_date

  data_type: 'integer'
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "ncbi_taxon_id",
  { data_type => "integer", is_nullable => 0 },
  "gene_id",
  { data_type => "integer", is_nullable => 0 },
  "discontinued_gene_id",
  { data_type => "integer", is_nullable => 1 },
  "discontinued_symbol",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "discontinue_date",
  { data_type => "integer", is_nullable => 1 },
);


# Created by DBIx::Class::Schema::Loader v0.07035 @ 2013-07-18 14:44:41
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:it1UDtKpEpdGmLtn47mweA


# You can replace this text with custom code or comments, and it will be preserved on regeneration
__PACKAGE__->meta->make_immutable;
1;
