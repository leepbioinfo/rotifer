use utf8;
package Rotifer::DBIC::AnnotationDB::Result::Ontology;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Rotifer::DBIC::AnnotationDB::Result::Ontology

=cut

use strict;
use warnings;

use Moose;
use MooseX::NonMoose;
use MooseX::MarkAsMethods autoclean => 1;
extends 'DBIx::Class::Core';

=head1 TABLE: C<ontology>

=cut

__PACKAGE__->table("ontology");

=head1 ACCESSORS

=head2 ontology_id

  data_type: 'integer'
  is_auto_increment: 1
  is_nullable: 0

=head2 name

  data_type: 'varchar'
  is_nullable: 0
  size: 255

=head2 definition

  data_type: 'text'
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "ontology_id",
  { data_type => "integer", is_auto_increment => 1, is_nullable => 0 },
  "name",
  { data_type => "varchar", is_nullable => 0, size => 255 },
  "definition",
  { data_type => "text", is_nullable => 1 },
);

=head1 PRIMARY KEY

=over 4

=item * L</ontology_id>

=back

=cut

__PACKAGE__->set_primary_key("ontology_id");

=head1 UNIQUE CONSTRAINTS

=head2 C<ontology_name_uniq>

=over 4

=item * L</name>

=back

=cut

__PACKAGE__->add_unique_constraint("ontology_name_uniq", ["name"]);

=head1 RELATIONS

=head2 term_paths

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::TermPath>

=cut

__PACKAGE__->has_many(
  "term_paths",
  "Rotifer::DBIC::AnnotationDB::Result::TermPath",
  { "foreign.ontology_id" => "self.ontology_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 term_relationships

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::TermRelationship>

=cut

__PACKAGE__->has_many(
  "term_relationships",
  "Rotifer::DBIC::AnnotationDB::Result::TermRelationship",
  { "foreign.ontology_id" => "self.ontology_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 terms

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::Term>

=cut

__PACKAGE__->has_many(
  "terms",
  "Rotifer::DBIC::AnnotationDB::Result::Term",
  { "foreign.ontology_id" => "self.ontology_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);


# Created by DBIx::Class::Schema::Loader v0.07033 @ 2013-05-25 23:07:20
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:MY9wUn/bZrybO9w8Hn7/Jg


# You can replace this text with custom code or comments, and it will be preserved on regeneration
__PACKAGE__->meta->make_immutable;
1;
