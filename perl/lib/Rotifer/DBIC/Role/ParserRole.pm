# POD documentation - main docs before the code

=head1 NAME

Rotifer::DBIC::Role::ParserRole - basic attributes, methods of a parser 

=head1 SYNOPSIS

  # Implementing a new parser

  with 'Rotifer::DBIC::Role::ParserRole';

  sub load {
    ...do something...
  }

=head1 DESCRIPTION

This moudle provides basic method iplementations and attributes
and defines an interface for Rotifer::DBIC parsers.

Use it as a base when implementing your own parser.

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::DBIC::Role::ParserRole;

use strict;
use warnings;
#use autodie qw(:all);
use Carp::Clan qr/^Rotifer::DBIC/;
use Moose::Role;
with 'Rotifer::DBIC::Role::BasePluginRole';

=head2 ATTRIBUTES / ACCESSORS

=head2 batch_size

 Usage   : Rotifer::DBIC::ExampleDB::Parser->new(batch_size => 6000)
 Function: number of sequences to process per batch
 Value   : integer
 Default :

=cut

has 'batch_size'    => (
    is      => "rw",
    isa     => "Int",
    default => 5000
    );

=head2 METHODS

=head2 load

 Title   : load
 Usage   : @biodata = $parser->load(@files)
 Function: process and load data from all input files
 Returns : array of Rotifer::DBIC::ExampleDB::Result::Biodata
 Args    : list of file names and/or file handlers

=cut

requires 'load';

1;
