# POD documentation - main docs before the code

=head1 NAME

Rotifer::DBIC::Role::BasePluginRole - 
 general properties/methods for Rotifer::DBIC components 

=head1 SYNOPSIS

  # Implementing a new parser

  with 'Rotifer::DBIC::Role::BasePluginRole';

  sub load {
    ...do something...
  }

=head1 DESCRIPTION

Parsers, writers and similar plugins are some of the ppossible
components of a Rotifer::DBIC schema. This module provides the
basic method iplementations, attributes and also defines the
methods a Rotifer::DBIC plugin should support.

Use it as a base when implementing your own plugins.

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::DBIC::Role::BasePluginRole;

use strict;
use warnings;
#use autodie qw(:all);
use Carp::Clan qr/^Rotifer::DBIC/;
use Moose::Role;

=head2 ATTRIBUTES / ACCESSORS

=head2 schema

 Usage   : Rotifer::DBIC::ExampleDB::Writer->new(schema => $rdbic)
 Function: set database connection
 Value   : Rotifer::DBIC::AnnotationDB schema object

=cut

has 'schema'   => (
    is       => "rw",
    does     => "Rotifer::DBIC::Role::BaseSchemaRole",
    required => 1
    );

=head2 isa_transaction_manager

 Usage   : $writer->isa_transaction_manager(1)
 Function: true if the plugin will take care of all its transactions
 Value   : boolean
 Default : false

=cut

has 'isa_transaction_manager' => (
    is  => "ro",
    isa => 'Bool',
    default => 0,
    );

=head2 METHODS

No methods are required by this role.

=cut

1;
