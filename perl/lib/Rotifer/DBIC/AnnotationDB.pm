use utf8;
package Rotifer::DBIC::AnnotationDB;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

use Moose;
use MooseX::MarkAsMethods autoclean => 1;
extends 'DBIx::Class::Schema';

__PACKAGE__->load_namespaces;


# Created by DBIx::Class::Schema::Loader v0.07033 @ 2013-05-02 15:34:59
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:mOGiearsFU/2WRWFL5hb3w

# POD documentation - main docs before the code

=head1 NAME

Rotifer::DBIC::AnnotationDB - Rotifer database for sequence annotation

=head1 SYNOPSIS

  # Using this library through rannotationDB
  rannotationDB -if gi -of phyloprofile input.gi

  # If you would like to program:
  #
  # Creating a new parser
  use Rotifer::DBIC::AnnotationDB;
  my $db = Rotifer::DBIC::AnnotationDB->new(%opts);
  my @biodata = $db->resultset("Biodata")->search({ });

=head1 DESCRIPTION

Rotifer::DBIC::::AnnotationDB is a SQL database backend for storing
sequence annotations. Inspired by BioSQL and Chado, this database
stores data in a set of tables that allow for flexibile sequence
annotation and vocabulary.

Bundled with the database schema, a set of parsers and writers allow
loading and querying the database from command line in an easy way.

The rannotationDB command-line application implements an interface to
Rotifer::DBIC databases that supports Rotifer::DBIC::AnnotationDB
and can be used as a prototype for any application that needs to use
the backend. 

=head2 Data Types in Rotifer::DBIC

For the purpose of correctly modelling your data while implementing new parsers,
it is useful to understand how Rotifer::DBIC represents data internally.

Essentially, there are three separate data environments in Rotifer::DBIC:

=over

=item (1) Data entries

Data loaded by the user using a Rotifer::DBIC::AnnotationDB parser 
other than the ontology parsers, including sequences, external 
database references, user annotations, etc.

=item (2) Datasets

Datasets are collections of data entries or other datasets.

Note that every time data is loaded, AnnotationDB ensures a new 
dataset is created to group the new data.

Projects are a special type of dataset, supported by rannotationDB,
that aggregates datasets and may be used in searches. 

=item (3) Ontologies

Just like Chado or BioSQL, AnnotationDB uses ontologies as the source
for a vocabulary of terms used to describe and type data and its 
annotations tags. For example, if one would like to associate a 
comment to some datum (e.g. a sequence), he/she will have to first
add the word 'comment' as a keyword in the current ontology.

Note to developers: parsers are expected to create ontologies 
automatically when previously loaded ontologies are not explicitly
chosen by the user. 

=back

=head1 DEPENDENCIES

Consumed roles, inherited classes and imported methods.

=over

=item Rotifer::DBIC::Role::BaseSchemaRole

=item Rotifer::DBIC::AnnotationDB::Role::HasDatasetsRole

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

# Dependencies
use Carp::Clan qr/^Rotifer::DBIC::AnnotationDB/;
use Scalar::Util qw(blessed);
use Rotifer;
with qw( 
 Rotifer::DBIC::Role::BaseSchemaRole
);
our $VERSION = $Rotifer::VERSION;

=head2 ATTRIBUTES / ACCESSORS

=head2 datasets

 Usage   : $parser->datasets([ $set1, $set2 ])
 Function: set/select datasets
 Value   : arrayref of Datasets
 Default : reference to empty array

=cut

has 'datasets' => (
    is      => "rw", 
    isa     => "ArrayRef[Object]",
    lazy    => 1,
    default => sub { [] },
    );

=head2 METHODS

See documentation for DBIc::Class::Schema and Rotifer::DBIC::Role::BaseSchemaRole.

=cut

sub deploy {
    my $self = shift;
    my @ret = $self->next::method(@_);

    # Ugly hack for BioSQL compatibility on PostgreSQL
    #  This will be necessary until I learn how to force
    #  DBIx::Class::Schema::Loader to transfer sequence names
    #  or t move away from BioSQL (best option!)
    if ($self->storage->sqlt_type eq 'PostgreSQL') {
	$self->storage->debug($self->debug > 1 ? 1 : 0);
	my @stuff = $self->storage->dbh_do(
	    sub {
		my ($storage, $dbh, @args) = @_;
		my %fixed = ();
		my %biosql = map { ($_,1) } qw(dbxref dbxref_qualifier_value ontology term term_dbxref term_path term_relationship term_synonym);
		foreach my $source_name ($self->sources) {
		    my $source = $self->source($source_name);
		    my $name = $source->name;
		    next if (exists $fixed{$name}); # Why does this loop process a source more than once?
		    my @auto = grep { exists $source->column_info($_)->{is_auto_increment} } $source->columns;
		    if (scalar @auto) {
			my $q = "ALTER SEQUENCE ${name}_$auto[0]_seq RENAME TO ${name}_pk_seq";
			print STDERR "$q\n" if ($self->debug);
			$dbh->do($q);
		    }
		    $fixed{$name}++;
		}
	    });
    } # if ($self->storage->sqlt_type eq 'PostgreSQL')

    return @ret;
}

sub deployed {
    return scalar(grep { /biodata/i } shift->storage->dbh->tables) ? 1 : 0;
}

__PACKAGE__->meta->make_immutable(inline_constructor => 0);
1;
