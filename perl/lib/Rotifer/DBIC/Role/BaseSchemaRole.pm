# POD documentation - main docs before the code

=head1 NAME

Rotifer::DBIC::Role::BaseSchemaRole -
 basic attributes and methods for a Rotifer::DBIC schema 

=head1 SYNOPSIS

  # Implementing a new database

  with 'Rotifer::DBIC::Role::BaseSchemaRole';

  sub load {
    ...do something...
  }

  sub write {
    ...do something...
  }

=head1 DESCRIPTION

This moudle defines basic method iplementations or requirements and
basic attributes of an interface for Rotifer::DBIC databases.

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

package Rotifer::DBIC::Role::BaseSchemaRole;

use strict;
use warnings;
#use autodie qw(:all);
use Carp::Clan qr/^Rotifer::DBIC/;
use Scalar::Util qw(blessed);
use Moose::Role;

=head2 ATTRIBUTES / ACCESSORS

=head2 debug

 Usage   : Rotifer::DBIC::AnnotationDB::Writer->new(debug => 1)
 Function: set debug messages level
 Value   : integer

=cut

has 'debug'    => (
    is      => "rw",
    isa     => "Str",
    default => 0
    );

=head2 METHODS

=head2 deployed

 Title   : deployed
 Usage   : my $ok = $schema->deployed
 Function: check whether an instance of this schema exists
 Returns : boolean
 Args    : none

=cut

requires 'deployed';

=head2 write

 Title   : write
 Usage   : @biodata = $parser->write("fasta",\%opts, @ids)
 Function: query database and dump data to some external format
 Returns : boolean, true on success
 Args    : output format, writer options and a list of file names,
           file handlers or whatever input is supported by the writer

=cut

sub write {
    my ($self, $format, $opts, @args) = @_;
    my $worker_class = blessed($self) ? ref($self) : $self;
    $worker_class .= "::Writer";
    eval "require $worker_class" or die $@;
    $opts->{'schema'}  = $self;

    # Call writer plugin
    die "Rotifer::DBIC writer factories should be based on MooseX::AbstractFactory"
	unless ($worker_class->can("create"));
    my $worker = $worker_class->create($format, $opts);
    die "Writers should consume Rotifer::DBIC::Role::WriterRole" unless $worker->does("Rotifer::DBIC::Role::WriterRole");
    my $guard = $self->txn_scope_guard unless ($worker->isa_transaction_manager);
    my $status = $worker->write(@args);
    $guard->commit if (defined $guard);

    return $status;
}

=head2 load

 Title   : load
 Usage   : @biodata = $parser->load("gi", \%opts, @files)
 Function: process and load data into database
 Returns : array of Rotifer::DBIC::AnnotationDB::Result::Biodata
 Args    : input format, parser options and a list of file names,
           file handlers or whatever is supported by the parser

=cut

sub load {
    my ($self, $format, $opts, @args) = @_;
    my $worker_class = blessed($self) ? ref($self) : $self;
    $worker_class .= "::Parser";
    eval "require $worker_class" or die $@;
    $opts->{'schema'}  = $self;

    # Create a new instance of Rotifer::DBIC
    if (!$self->deployed) {
	$self->storage->txn_do(
	    sub {
		$self->deploy({},$ENV{"PWD"});
	    });
    }

    # Call parser plugin
    die "Rotifer::DBIC parser factories should be based on MooseX::AbstractFactory"
	unless ($worker_class->can("create"));
    my $worker = $worker_class->create($format, $opts);
    die "Parsers should consume Rotifer::DBIC::Role::ParserRole" unless $worker->does("Rotifer::DBIC::Role::ParserRole");
    my $guard = $self->txn_scope_guard unless ($worker->isa_transaction_manager);
    my @val = $worker->load(@args);
    $guard->commit if (defined $guard);

    return @val;
}

# Make perl happy
1;
