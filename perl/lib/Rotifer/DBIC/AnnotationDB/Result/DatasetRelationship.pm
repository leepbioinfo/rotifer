use utf8;
package Rotifer::DBIC::AnnotationDB::Result::DatasetRelationship;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Rotifer::DBIC::AnnotationDB::Result::DatasetRelationship

=cut

use strict;
use warnings;

use Moose;
use MooseX::NonMoose;
use MooseX::MarkAsMethods autoclean => 1;
extends 'DBIx::Class::Core';

=head1 TABLE: C<dataset_relationship>

=cut

__PACKAGE__->table("dataset_relationship");

=head1 ACCESSORS

=head2 dataset_relationship_id

  data_type: 'integer'
  is_nullable: 0

=head2 subject_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 1

=head2 object_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "dataset_relationship_id",
  { data_type => "integer", is_nullable => 0 },
  "subject_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 1 },
  "object_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 1 },
);

=head1 PRIMARY KEY

=over 4

=item * L</dataset_relationship_id>

=back

=cut

__PACKAGE__->set_primary_key("dataset_relationship_id");

=head1 UNIQUE CONSTRAINTS

=head2 C<dataset_relationship_uniq>

=over 4

=item * L</subject_id>

=item * L</object_id>

=back

=cut

__PACKAGE__->add_unique_constraint("dataset_relationship_uniq", ["subject_id", "object_id"]);

=head1 RELATIONS

=head2 object

Type: belongs_to

Related object: L<Rotifer::DBIC::AnnotationDB::Result::Dataset>

=cut

__PACKAGE__->belongs_to(
  "object",
  "Rotifer::DBIC::AnnotationDB::Result::Dataset",
  { dataset_id => "object_id" },
  {
    is_deferrable => 1,
    join_type     => "LEFT",
    on_delete     => "CASCADE",
    on_update     => "CASCADE",
  },
);

=head2 subject

Type: belongs_to

Related object: L<Rotifer::DBIC::AnnotationDB::Result::Dataset>

=cut

__PACKAGE__->belongs_to(
  "subject",
  "Rotifer::DBIC::AnnotationDB::Result::Dataset",
  { dataset_id => "subject_id" },
  {
    is_deferrable => 1,
    join_type     => "LEFT",
    on_delete     => "CASCADE",
    on_update     => "CASCADE",
  },
);


# Created by DBIx::Class::Schema::Loader v0.07033 @ 2013-06-17 01:39:13
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:lgmifNOF00sehq9cDiSr6g


# You can replace this text with custom code or comments, and it will be preserved on regeneration
__PACKAGE__->meta->make_immutable;
1;
