# POD documentation - main docs before the code

=head1 NAME

Rotifer::DBIC::AnnotationDB::Role::SearchRole - 
  basic attributes and methods of an AnnotationDB search module 

=head1 SYNOPSIS

  # Implementing a new writer

  with 'Rotifer::DBIC::AnnotationDB::Role::SearchRole';

=head1 DESCRIPTION

This moudle provides basic method iplementations and attributes
and defines an interface for Rotifer::DBIC::AnnotationDB search
modules.

Use it as a base when implementing your own writer.

=head1 DEPENDENCIES

Consumed roles, inherited classes and imported methods.

=over

=item Rotifer::DBIC::Role::SearchRole

=item Scalar::Util

=back

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::DBIC::AnnotationDB::Role::SearchRole;

use strict;
use warnings;
#use autodie qw(:all);
use Carp::Clan qr/^Rotifer::DBIC::AnnotationDB/;
use Moose::Role;
with qw(Rotifer::DBIC::Role::SearchRole);

=head2 ATTRIBUTES / ACCESSORS

See dependencies (above) for attributes, etc.

=head2 METHODS

See dependencies (above) for attributes, etc.

=cut

1;
