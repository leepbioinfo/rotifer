# POD documentation - main docs before the code

=head1 NAME

Rotifer::DBIC::AnnotationDB::Role::ParserRole - 
  basic attributes, methods for AnnotationDB parsers 

=head1 SYNOPSIS

  # Implementing a new parser

  with 'Rotifer::DBIC::AnnotationDB::Role::ParserRole';

=head1 DESCRIPTION

This moudle provides basic method iplementations and attributes
and defines an interface for Rotifer::DBIC::AnnotationDB parsers.

Use it as a base when implementing your own parser.

=head1 DEPENDENCIES

Consumed roles, inherited classes and imported methods.

=over

=item Rotifer::DBIC::Role::ParserRole

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

package Rotifer::DBIC::AnnotationDB::Role::ParserRole;

use strict;
use warnings;
#use autodie qw(:all);
use Carp::Clan qr/^Rotifer::DBIC::AnnotationDB/;
use Moose::Role;
with 'Rotifer::DBIC::Role::ParserRole';

=head2 ATTRIBUTES / ACCESSORS

=cut

=head2 parser_ontology

 Usage   : $parser->parser_ontology
 Function: ontology for parser-specific annotations tags
 Returns : array reference
 Builder : _default_parser_ontology
 Trigger : _parser_ontology_trigger

=cut

has 'parser_ontology' => (
    is       => "rw", 
    isa      => "Str|Object",
    lazy     => 1,
    builder  => '_default_parser_ontology',
    trigger  => \&_parser_ontology_trigger,
    );

sub _default_parser_ontology {
    return shift->schema->resultset('Ontology')->
	find_or_create({ name => "AnnotationDB Parser Tags" });
}

sub _parser_ontology_trigger {
    my ($self, $new, $old) = @_;
    unless (blessed $new && $new->isa("Rotifer::DBIC::AnnotationDB::Result::Ontology")) {
	$new = { name => $new } unless UNIVERSAL::isa($new,"HASH");
	my $onto = $self->schema->resultset('Ontology')->update_or_create($new);
	if (defined $onto && blessed $onto && $onto->isa("Rotifer::DBIC::AnnotationDB::Result::Ontology")) {
	    $self->parser_ontology($onto);
	} else {
	    my $name = ref $new ? $new->{name} : $new;
	    croak "Failed to load ontology $name: $@";
	} 
    }
}

=head2 tags

 Usage   : $parser->tags
 Function: annotation tags used by this parser
 Returns : Rotifer::DBIC::AnnotationDB::Result::Term
 Builder : _default_tags

 The following tags will be added to the default vocabulary
 unless present in one of the other ontologies listed by
 $parser->parser_ontology (see ATTRIBUTES above):

 fasta   : regular file containing FASTA sequences

=cut

has 'tags' => (
    is       => "ro", 
    isa      => "HashRef[Object]",
    lazy     => 1,
    builder  => '_default_tags',
    init_arg => undef,
    );

sub _default_tags {
    my $self = shift;
    return {};
}

sub _hash2tags {
    my ($self, $onto, @tags) = @_;
    my $trs = $self->schema->resultset("Term");
    my $tsrs = $self->schema->resultset('TermSynonym');

    # Turning tags into objects
    my %tags = ();
    foreach my $hash (@tags) {
	$hash->{is_obsolete} = 0; # Make sure no obsolete tags are used
	my $synonyms = exists $hash->{synonyms} ? delete $hash->{synonyms} : ();
	$hash->{ontology} = $onto;
	my $term = $trs->find_or_create($hash);
	map { $term->find_or_create_related('term_synonyms',{ synonym => $_ }) } @$synonyms;
	$tags{$hash->{name}} = $term;
    }

    return \%tags;
}

=head2 METHODS

=head2 load

 Title   : load
 Usage   : @biodata = $parser->load(@files)
 Function: process and load data from all input files
 Returns : array of Rotifer::DBIC::ExampleDB::Result::Biodata
 Args    : list of file names and/or file handlers

 The method overrides RotifeR::DBIC::Role::ParserRole
 to make sure each load() call creates a new dataset
 and adds this dataset to the schema.

=cut

around 'load' => sub {
    my ($orig, $self, @args) = @_;

    # Add new dataset to schema
    my $nset = undef;
    if (!scalar @{ $self->schema->datasets }) {
	$nset = $self->_auto_create_dataset;
	$self->schema->datasets([ $nset ]);
    }

    # Run real plugin method
    my @data = $self->$orig(@args);

    # Remove dataset if nothing was added to it
    if (defined $nset) {
	my $keep = 0;
	foreach my $relname ($nset->result_source->relationships) {
	    next unless ($relname =~ /_datasets$/);
	    my $count = $nset->search_related($relname)->count;
	    $keep = 1 if ($count);
	}
	if (!$keep) {
	    $self->schema->datasets([]);
	    $nset->delete;
	}
    }

    return @data;
};

# Helper function
sub _auto_create_dataset {
    my ($self) = @_;

    my $drs = $self->schema->resultset("Dataset");
    my $hash = {
	name        => Data::UUID->new->create_str,
	authority   => $ENV{USER},
	loaded_by   => ref $self,
	loaded_from => $ENV{PWD},
    };
    my $dataset = $drs->find_or_create($hash,{ key => "dataset_uniq" });

    return $dataset;
}

1;
