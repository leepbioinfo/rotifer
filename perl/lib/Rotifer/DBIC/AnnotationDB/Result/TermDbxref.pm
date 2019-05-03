use utf8;
package Rotifer::DBIC::AnnotationDB::Result::TermDbxref;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Rotifer::DBIC::AnnotationDB::Result::TermDbxref

=cut

use strict;
use warnings;

use Moose;
use MooseX::NonMoose;
use MooseX::MarkAsMethods autoclean => 1;
extends 'DBIx::Class::Core';

=head1 TABLE: C<term_dbxref>

=cut

__PACKAGE__->table("term_dbxref");

=head1 ACCESSORS

=head2 term_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=head2 dbxref_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=head2 rank

  data_type: 'smallint'
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "term_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
  "dbxref_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
  "rank",
  { data_type => "smallint", is_nullable => 1 },
);

=head1 PRIMARY KEY

=over 4

=item * L</term_id>

=item * L</dbxref_id>

=back

=cut

__PACKAGE__->set_primary_key("term_id", "dbxref_id");

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
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:+NED8hbk91y0LugezTPvvw


# You can replace this text with custom code or comments, and it will be preserved on regeneration
__PACKAGE__->meta->make_immutable;
1;
