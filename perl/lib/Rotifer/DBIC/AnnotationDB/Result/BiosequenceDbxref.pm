use utf8;
package Rotifer::DBIC::AnnotationDB::Result::BiosequenceDbxref;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Rotifer::DBIC::AnnotationDB::Result::BiosequenceDbxref

=cut

use strict;
use warnings;

use Moose;
use MooseX::NonMoose;
use MooseX::MarkAsMethods autoclean => 1;
extends 'DBIx::Class::Core';

=head1 TABLE: C<biosequence_dbxref>

=cut

__PACKAGE__->table("biosequence_dbxref");

=head1 ACCESSORS

=head2 biosequence_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=head2 dbxref_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=cut

__PACKAGE__->add_columns(
  "biosequence_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
  "dbxref_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
);

=head1 PRIMARY KEY

=over 4

=item * L</biosequence_id>

=item * L</dbxref_id>

=back

=cut

__PACKAGE__->set_primary_key("biosequence_id", "dbxref_id");

=head1 RELATIONS

=head2 biosequence

Type: belongs_to

Related object: L<Rotifer::DBIC::AnnotationDB::Result::Biosequence>

=cut

__PACKAGE__->belongs_to(
  "biosequence",
  "Rotifer::DBIC::AnnotationDB::Result::Biosequence",
  { biosequence_id => "biosequence_id" },
  { is_deferrable => 1, on_delete => "CASCADE", on_update => "CASCADE" },
);

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


# Created by DBIx::Class::Schema::Loader v0.07033 @ 2013-05-07 16:31:56
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:aIR45rT92HS1PzmU1gFu+g


# You can replace this text with custom code or comments, and it will be preserved on regeneration
__PACKAGE__->meta->make_immutable;
1;
