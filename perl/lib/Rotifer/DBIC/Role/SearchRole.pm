# POD documentation - main docs before the code

=head1 NAME

Rotifer::DBIC::Role::WriterRole - attributes and methods of a Rotifer::DBIC writer 

=head1 SYNOPSIS

  # Implementing a new parser

  with 'Rotifer::DBIC::Role::WriterRole';

  sub load {
    ...do something...
  }

=head1 DESCRIPTION

This moudle provides basic method iplementations and attributes
and defines an interface for Rotifer::DBIC::AnnotationDB parsers.

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

package Rotifer::DBIC::Role::WriterRole;

use strict;
use warnings;
#use autodie qw(:all);
use Carp::Clan qr/^Rotifer::DBIC/;
use Moose::Role;
with 'Rotifer::DBIC::Role::BasePluginRole';

=head2 ATTRIBUTES / ACCESSORS

=head2 input_format

 Usage   : $db->input_format("gi")
 Function: type of queries given as input for searches
 Value   : string
 Default : undefined

 Acceptable values for this attribute are implementation
 specific.

=cut

has 'input_format' => (
    is      => "rw",
    isa     => "Str",
    );

=head2 supported_input_formats

 Usage   : $writer->supported_input_formats("fasta");
 Function: list of input formats that might be used or
           from which we can extract queries to search for
 Returns : hash reference
 Value   : none, attribute is read only

 Writers should are expected to override the default value
 for this attribute

=cut

has 'supported_input_formats' => (
    is      => 'ro',
    isa     => 'HashRef',
    default => sub { {} },
    );

=head2 METHODS

=head2 supports

 Title   : supports
 Usage   : @biodata = $parser->supports($format)
 Function: test whether the writer can extract queries
           from some input format
 Returns : boolean
 Args    : string

=cut

sub supports {
    my ($self,$format) = @_;
    my @keys = keys %{ $self->supported_input_formats };
    die "Writer module did not override default value for attribute 'supported_input_formats'"
	unless (scalar @keys);
    return 1 if (exists $self->supported_input_formats->{$format});
    my @t = grep { $format =~ /^${_}::\S+$/ } @keys;
    return scalar @t ? 1 : 0;
}

=head2 write

 Title   : write
 Usage   : @biodata = $parser->write(@files)
 Function: process and write data from all input files
 Returns : array of Rotifer::DBIC::AnnotationDB::Result::Biodata
 Args    : list of ids, file names, file handlers or objects

=cut

requires qw(write);

before "write" => sub {
    my $self = shift;
    if (scalar @_) {
	my $format = $self->input_format;
	if (!defined $format) {
	    die "To use input data as queries you must name an input format.";
	}
	my $name = blessed $self ? ref($self) : $self;
	die "Writer $name doesn't know how to process input format $format"
	    unless $self->supports($format);
    }
};

1;
