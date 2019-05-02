use utf8;
package Rotifer::DBIC::AnnotationDB::Result::TermPath;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Rotifer::DBIC::AnnotationDB::Result::TermPath

=cut

use strict;
use warnings;

use Moose;
use MooseX::NonMoose;
use MooseX::MarkAsMethods autoclean => 1;
extends 'DBIx::Class::Core';

=head1 TABLE: C<term_path>

=cut

__PACKAGE__->table("term_path");

=head1 ACCESSORS

=head2 term_path_id

  data_type: 'integer'
  is_auto_increment: 1
  is_nullable: 0

=head2 subject_term_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=head2 predicate_term_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=head2 object_term_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=head2 ontology_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=head2 distance

  data_type: 'integer'
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "term_path_id",
  { data_type => "integer", is_auto_increment => 1, is_nullable => 0 },
  "subject_term_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
  "predicate_term_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
  "object_term_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
  "ontology_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
  "distance",
  { data_type => "integer", is_nullable => 1 },
);

=head1 PRIMARY KEY

=over 4

=item * L</term_path_id>

=back

=cut

__PACKAGE__->set_primary_key("term_path_id");

=head1 UNIQUE CONSTRAINTS

=head2 C<term_path_subject_term_id_key>

=over 4

=item * L</subject_term_id>

=item * L</predicate_term_id>

=item * L</object_term_id>

=item * L</ontology_id>

=item * L</distance>

=back

=cut

__PACKAGE__->add_unique_constraint(
  "term_path_subject_term_id_key",
  [
    "subject_term_id",
    "predicate_term_id",
    "object_term_id",
    "ontology_id",
    "distance",
  ],
);

=head1 RELATIONS

=head2 object_term

Type: belongs_to

Related object: L<Rotifer::DBIC::AnnotationDB::Result::Term>

=cut

__PACKAGE__->belongs_to(
  "object_term",
  "Rotifer::DBIC::AnnotationDB::Result::Term",
  { term_id => "object_term_id" },
  { is_deferrable => 1, on_delete => "CASCADE", on_update => "CASCADE" },
);

=head2 ontology

Type: belongs_to

Related object: L<Rotifer::DBIC::AnnotationDB::Result::Ontology>

=cut

__PACKAGE__->belongs_to(
  "ontology",
  "Rotifer::DBIC::AnnotationDB::Result::Ontology",
  { ontology_id => "ontology_id" },
  { is_deferrable => 1, on_delete => "CASCADE", on_update => "CASCADE" },
);

=head2 predicate_term

Type: belongs_to

Related object: L<Rotifer::DBIC::AnnotationDB::Result::Term>

=cut

__PACKAGE__->belongs_to(
  "predicate_term",
  "Rotifer::DBIC::AnnotationDB::Result::Term",
  { term_id => "predicate_term_id" },
  { is_deferrable => 1, on_delete => "CASCADE", on_update => "CASCADE" },
);

=head2 subject_term

Type: belongs_to

Related object: L<Rotifer::DBIC::AnnotationDB::Result::Term>

=cut

__PACKAGE__->belongs_to(
  "subject_term",
  "Rotifer::DBIC::AnnotationDB::Result::Term",
  { term_id => "subject_term_id" },
  { is_deferrable => 1, on_delete => "CASCADE", on_update => "CASCADE" },
);


# Created by DBIx::Class::Schema::Loader v0.07033 @ 2013-06-07 04:40:13
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:+Gc/dHdVF3uvJ3YUWb/cug


# You can replace this text with custom code or comments, and it will be preserved on regeneration
__PACKAGE__->meta->make_immutable;
1;
