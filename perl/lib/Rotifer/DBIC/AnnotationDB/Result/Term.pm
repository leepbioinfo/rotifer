use utf8;
package Rotifer::DBIC::AnnotationDB::Result::Term;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Rotifer::DBIC::AnnotationDB::Result::Term

=cut

use strict;
use warnings;

use Moose;
use MooseX::NonMoose;
use MooseX::MarkAsMethods autoclean => 1;
extends 'DBIx::Class::Core';

=head1 TABLE: C<term>

=cut

__PACKAGE__->table("term");

=head1 ACCESSORS

=head2 term_id

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

=head2 identifier

  data_type: 'varchar'
  is_nullable: 1
  size: 40

=head2 is_obsolete

  data_type: 'char'
  default_value: 0
  is_nullable: 1
  size: 1

=head2 ontology_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=cut

__PACKAGE__->add_columns(
  "term_id",
  { data_type => "integer", is_auto_increment => 1, is_nullable => 0 },
  "name",
  { data_type => "varchar", is_nullable => 0, size => 255 },
  "definition",
  { data_type => "text", is_nullable => 1 },
  "identifier",
  { data_type => "varchar", is_nullable => 1, size => 40 },
  "is_obsolete",
  { data_type => "char", default_value => 0, is_nullable => 1, size => 1 },
  "ontology_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
);

=head1 PRIMARY KEY

=over 4

=item * L</term_id>

=back

=cut

__PACKAGE__->set_primary_key("term_id");

=head1 UNIQUE CONSTRAINTS

=head2 C<term_identifier_uniq>

=over 4

=item * L</identifier>

=back

=cut

__PACKAGE__->add_unique_constraint("term_identifier_uniq", ["identifier"]);

=head2 C<term_uniq>

=over 4

=item * L</name>

=item * L</ontology_id>

=item * L</is_obsolete>

=back

=cut

__PACKAGE__->add_unique_constraint("term_uniq", ["name", "ontology_id", "is_obsolete"]);

=head1 RELATIONS

=head2 biodata_attribute_values

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::BiodataAttributeValue>

=cut

__PACKAGE__->has_many(
  "biodata_attribute_values",
  "Rotifer::DBIC::AnnotationDB::Result::BiodataAttributeValue",
  { "foreign.term_id" => "self.term_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 biodata_relationship_attribute_values

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::BiodataRelationshipAttributeValue>

=cut

__PACKAGE__->has_many(
  "biodata_relationship_attribute_values",
  "Rotifer::DBIC::AnnotationDB::Result::BiodataRelationshipAttributeValue",
  { "foreign.term_id" => "self.term_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 biodata_relationships

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::BiodataRelationship>

=cut

__PACKAGE__->has_many(
  "biodata_relationships",
  "Rotifer::DBIC::AnnotationDB::Result::BiodataRelationship",
  { "foreign.predicate_id" => "self.term_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 biodatas

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::Biodata>

=cut

__PACKAGE__->has_many(
  "biodatas",
  "Rotifer::DBIC::AnnotationDB::Result::Biodata",
  { "foreign.term_id" => "self.term_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 biosequences

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::Biosequence>

=cut

__PACKAGE__->has_many(
  "biosequences",
  "Rotifer::DBIC::AnnotationDB::Result::Biosequence",
  { "foreign.term_id" => "self.term_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 dataset_attribute_values

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::DatasetAttributeValue>

=cut

__PACKAGE__->has_many(
  "dataset_attribute_values",
  "Rotifer::DBIC::AnnotationDB::Result::DatasetAttributeValue",
  { "foreign.term_id" => "self.term_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 dbxref_qualifier_values

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::DbxrefQualifierValue>

=cut

__PACKAGE__->has_many(
  "dbxref_qualifier_values",
  "Rotifer::DBIC::AnnotationDB::Result::DbxrefQualifierValue",
  { "foreign.term_id" => "self.term_id" },
  { cascade_copy => 0, cascade_delete => 0 },
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

=head2 organism_attribute_values

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::OrganismAttributeValue>

=cut

__PACKAGE__->has_many(
  "organism_attribute_values",
  "Rotifer::DBIC::AnnotationDB::Result::OrganismAttributeValue",
  { "foreign.term_id" => "self.term_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 term_dbxrefs

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::TermDbxref>

=cut

__PACKAGE__->has_many(
  "term_dbxrefs",
  "Rotifer::DBIC::AnnotationDB::Result::TermDbxref",
  { "foreign.term_id" => "self.term_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 term_path_object_terms

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::TermPath>

=cut

__PACKAGE__->has_many(
  "term_path_object_terms",
  "Rotifer::DBIC::AnnotationDB::Result::TermPath",
  { "foreign.object_term_id" => "self.term_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 term_path_predicate_terms

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::TermPath>

=cut

__PACKAGE__->has_many(
  "term_path_predicate_terms",
  "Rotifer::DBIC::AnnotationDB::Result::TermPath",
  { "foreign.predicate_term_id" => "self.term_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 term_path_subject_terms

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::TermPath>

=cut

__PACKAGE__->has_many(
  "term_path_subject_terms",
  "Rotifer::DBIC::AnnotationDB::Result::TermPath",
  { "foreign.subject_term_id" => "self.term_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 term_relationship_object_terms

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::TermRelationship>

=cut

__PACKAGE__->has_many(
  "term_relationship_object_terms",
  "Rotifer::DBIC::AnnotationDB::Result::TermRelationship",
  { "foreign.object_term_id" => "self.term_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 term_relationship_predicate_terms

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::TermRelationship>

=cut

__PACKAGE__->has_many(
  "term_relationship_predicate_terms",
  "Rotifer::DBIC::AnnotationDB::Result::TermRelationship",
  { "foreign.predicate_term_id" => "self.term_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 term_relationship_subject_terms

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::TermRelationship>

=cut

__PACKAGE__->has_many(
  "term_relationship_subject_terms",
  "Rotifer::DBIC::AnnotationDB::Result::TermRelationship",
  { "foreign.subject_term_id" => "self.term_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 term_synonyms

Type: has_many

Related object: L<Rotifer::DBIC::AnnotationDB::Result::TermSynonym>

=cut

__PACKAGE__->has_many(
  "term_synonyms",
  "Rotifer::DBIC::AnnotationDB::Result::TermSynonym",
  { "foreign.term_id" => "self.term_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);


# Created by DBIx::Class::Schema::Loader v0.07033 @ 2013-06-07 05:19:54
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:jVd8e+NME2+YCxakuhzV9A


# You can replace this text with custom code or comments, and it will be preserved on regeneration
__PACKAGE__->meta->make_immutable;
1;
