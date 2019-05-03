use utf8;
package Rotifer::DBIC::AnnotationDB::Result::OrganismDataset;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Rotifer::DBIC::AnnotationDB::Result::OrganismDataset

=cut

use strict;
use warnings;

use Moose;
use MooseX::NonMoose;
use MooseX::MarkAsMethods autoclean => 1;
extends 'DBIx::Class::Core';

=head1 TABLE: C<organism_dataset>

=cut

__PACKAGE__->table("organism_dataset");

=head1 ACCESSORS

=head2 organism_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=head2 dataset_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=cut

__PACKAGE__->add_columns(
  "organism_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
  "dataset_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
);

=head1 UNIQUE CONSTRAINTS

=head2 C<organism_dataset_uniq>

=over 4

=item * L</organism_id>

=item * L</dataset_id>

=back

=cut

__PACKAGE__->add_unique_constraint("organism_dataset_uniq", ["organism_id", "dataset_id"]);

=head1 RELATIONS

=head2 dataset

Type: belongs_to

Related object: L<Rotifer::DBIC::AnnotationDB::Result::Dataset>

=cut

__PACKAGE__->belongs_to(
  "dataset",
  "Rotifer::DBIC::AnnotationDB::Result::Dataset",
  { dataset_id => "dataset_id" },
  { is_deferrable => 1, on_delete => "CASCADE", on_update => "CASCADE" },
);

=head2 organism

Type: belongs_to

Related object: L<Rotifer::DBIC::AnnotationDB::Result::Organism>

=cut

__PACKAGE__->belongs_to(
  "organism",
  "Rotifer::DBIC::AnnotationDB::Result::Organism",
  { organism_id => "organism_id" },
  { is_deferrable => 1, on_delete => "CASCADE", on_update => "CASCADE" },
);


# Created by DBIx::Class::Schema::Loader v0.07033 @ 2013-06-17 13:33:01
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:33FxHd/JwdlXIhGzy+O0Gg


# You can replace this text with custom code or comments, and it will be preserved on regeneration
__PACKAGE__->meta->make_immutable;
1;
