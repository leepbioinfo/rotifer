use utf8;
package Rotifer::DBIC::AnnotationDB::Result::DbxrefQualifierValue;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Rotifer::DBIC::AnnotationDB::Result::DbxrefQualifierValue

=cut

use strict;
use warnings;

use Moose;
use MooseX::NonMoose;
use MooseX::MarkAsMethods autoclean => 1;
extends 'DBIx::Class::Core';

=head1 TABLE: C<dbxref_qualifier_value>

=cut

__PACKAGE__->table("dbxref_qualifier_value");

=head1 ACCESSORS

=head2 dbxref_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=head2 term_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=head2 rank

  data_type: 'smallint'
  default_value: 0
  is_nullable: 0

=head2 value

  data_type: 'text'
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "dbxref_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
  "term_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
  "rank",
  { data_type => "smallint", default_value => 0, is_nullable => 0 },
  "value",
  { data_type => "text", is_nullable => 1 },
);

=head1 PRIMARY KEY

=over 4

=item * L</dbxref_id>

=item * L</term_id>

=item * L</rank>

=back

=cut

__PACKAGE__->set_primary_key("dbxref_id", "term_id", "rank");

=head1 RELATIONS

=head2 dbxref

Type: belongs_to

Related object: L<Rotifer::DBIC::AnnotationDB::Result::Dbxref>

=cut

__PACKAGE__->belongs_to(
  "dbxref",
  "Rotifer::DBIC::AnnotationDB::Result::Dbxref",
  { dbxref_id => "dbxref_id" },
  { is_deferrable => 1, on_delete => "CASCADE", on_update => "CASCADE" },
);

=head2 term

Type: belongs_to

Related object: L<Rotifer::DBIC::AnnotationDB::Result::Term>

=cut

__PACKAGE__->belongs_to(
  "term",
  "Rotifer::DBIC::AnnotationDB::Result::Term",
  { term_id => "term_id" },
  { is_deferrable => 1, on_delete => "CASCADE", on_update => "CASCADE" },
);


# Created by DBIx::Class::Schema::Loader v0.07033 @ 2013-05-07 16:31:56
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:XcebhPCxAtKoFz6Ygr3Icw


# You can replace this text with custom code or comments, and it will be preserved on regeneration
__PACKAGE__->meta->make_immutable;
1;
