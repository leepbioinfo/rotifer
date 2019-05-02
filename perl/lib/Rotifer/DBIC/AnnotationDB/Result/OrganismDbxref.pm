use utf8;
package Rotifer::DBIC::AnnotationDB::Result::OrganismDbxref;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Rotifer::DBIC::AnnotationDB::Result::OrganismDbxref

=cut

use strict;
use warnings;

use Moose;
use MooseX::NonMoose;
use MooseX::MarkAsMethods autoclean => 1;
extends 'DBIx::Class::Core';

=head1 TABLE: C<organism_dbxref>

=cut

__PACKAGE__->table("organism_dbxref");

=head1 ACCESSORS

=head2 organism_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=head2 dbxref_id

  data_type: 'integer'
  is_foreign_key: 1
  is_nullable: 0

=cut

__PACKAGE__->add_columns(
  "organism_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
  "dbxref_id",
  { data_type => "integer", is_foreign_key => 1, is_nullable => 0 },
);

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
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:Zd+XmBgzZ/v2qWVfTKPj0A


# You can replace this text with custom code or comments, and it will be preserved on regeneration
__PACKAGE__->meta->make_immutable;
1;
