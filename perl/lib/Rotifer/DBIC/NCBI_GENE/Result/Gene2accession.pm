use utf8;
package Rotifer::DBIC::NCBI_GENE::Result::Gene2accession;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Rotifer::DBIC::NCBI_GENE::Result::Gene2accession

=cut

use strict;
use warnings;

use Moose;
use MooseX::NonMoose;
use MooseX::MarkAsMethods autoclean => 1;
extends 'DBIx::Class::Core';

=head1 TABLE: C<gene2accession>

=cut

__PACKAGE__->table("gene2accession");

=head1 ACCESSORS

=head2 ncbi_taxon_id

  data_type: 'integer'
  is_nullable: 0

=head2 gene_id

  data_type: 'integer'
  is_nullable: 0

=head2 rank

  data_type: 'integer'
  default_value: 0
  is_nullable: 0

=head2 status

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 rna_accession

  data_type: 'varchar'
  is_nullable: 1
  size: 64

=head2 rna_gi

  data_type: 'integer'
  is_nullable: 1

=head2 protein_accession

  data_type: 'varchar'
  is_nullable: 1
  size: 64

=head2 protein_gi

  data_type: 'integer'
  is_nullable: 1

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

=head2 assembly

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 peptide_accession

  data_type: 'varchar'
  is_nullable: 1
  size: 64

=head2 peptide_gi

  data_type: 'integer'
  is_nullable: 1

=head2 symbol

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=cut

__PACKAGE__->add_columns(
  "ncbi_taxon_id",
  { data_type => "integer", is_nullable => 0 },
  "gene_id",
  { data_type => "integer", is_nullable => 0 },
  "rank",
  { data_type => "integer", default_value => 0, is_nullable => 0 },
  "status",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "rna_accession",
  { data_type => "varchar", is_nullable => 1, size => 64 },
  "rna_gi",
  { data_type => "integer", is_nullable => 1 },
  "protein_accession",
  { data_type => "varchar", is_nullable => 1, size => 64 },
  "protein_gi",
  { data_type => "integer", is_nullable => 1 },
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
  "assembly",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "peptide_accession",
  { data_type => "varchar", is_nullable => 1, size => 64 },
  "peptide_gi",
  { data_type => "integer", is_nullable => 1 },
  "symbol",
  { data_type => "varchar", is_nullable => 1, size => 255 },
);


# Created by DBIx::Class::Schema::Loader v0.07035 @ 2013-07-18 14:44:41
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:Ccb+Rtzf6YO1mt96UNOLdA


# You can replace this text with custom code or comments, and it will be preserved on regeneration
__PACKAGE__->meta->make_immutable;
1;
