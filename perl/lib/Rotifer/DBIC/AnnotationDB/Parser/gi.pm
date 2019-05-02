# POD documentation - main docs before the code

=head1 NAME

Rotifer::DBIC::AnnotationDB::Parser::gi - load sequences from GI/accession lists 

=head1 SYNOPSIS

  # Simplest usage: GIs or accession from the command line
  #
  # Note: this will automatically create a graph
  #       and add all sequences to it

  rannotationDB -a load -if gi 123456 YP_003577

  # Loading from a file 

  rannotationDB -a load -if gi gi.txt

  # Same rannotationDB thing from a pipe

  cat gi.txt | rannotationDB -a load -if gi

  # Creating a new parser script
  use Rotifer::DBIC::AnnotationDB::Parser;
  my $parser = Rotifer::DBIC::AnnotationDB::Parser->create("gi");
  $parser->load(@ARGV);

=head1 DESCRIPTION

Rotifer::DBIC::AnnotationDB::Parser::gi will parse the command line and any files
given as arguments to extract lists of GIs and/or accessions. 

It will then retrieve all sequences from local and/or remote
databases and load these sequences to a Rotifer::DBIC::AnnotationDB database. 

=head1 EXTENDS

 Rotifer::DBIC::AnnotationDB::Parser::fasta

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::DBIC::AnnotationDB::Parser::gi;

use strict;
use warnings;
use autodie qw(:all);
use Carp::Clan qr/^Rotifer::DBIC::AnnotationDB/;
use File::Monitored qw(ttempfile);
use IO::String;
use Moose;
use Rotifer::DB qw(id2fasta);
extends 'Rotifer::DBIC::AnnotationDB::Parser::fasta';

=head2 ATTRIBUTES / ACCESSORS

This section list attributes not defined by the fundamental AnnotationDB parser's role.

=head2 blastdb

 Usage   : Rotifer::DBIC::AnnotationDB::Parser->new(blastdb => "pdb")
 Function: name of the local BLAST database to use while retrieving sequences

           This parameter is used by id2fasta to speed up sequence retrieval
           and to allow for automatic loading of identical sequences

           You have to make sure fastacmd and/or blastdbcmd are able to use the
           database name you provide otherwise it will be ignored.

           When in doubt, provide the fully qualified path to the database.

 Value   : string
 Default : nr

=cut

has 'blastdb' => (
    is   => 'rw',
    isa  => 'Str',
    default => 'nr',
    );

=head2 redundant

 Usage   : Rotifer::DBIC::AnnotationDB::Parser->new(redundant => 1)
 Function: load all identical sequences from local BLAST
           databases using fastacmd
 Value   : boolean
 Default : 0

=cut

has 'redundant' => (
    is   => 'rw',
    isa  => 'Bool',
    default => 0,
    );

=head2 tags

 Usage   : $parser->tags
 Function: annotation tags used by this parser
 Returns : Rotifer::DBIC::AnnotationDB::Result::Term
 Builder : _default_tags

 The following tags will be added to the default vocabulary
 unless present in one of the preferred ontologies:

  gi file : a file with one sequence ID per row

=cut

sub _default_tags {
    my ($self) = shift;

    my @tags = 
	(
	 {
	     name       => 'gi file',
	     definition => "Text file with a list of sequence identifiers (GI or accession numbers), one per row.",
	 },
	);

    return $self->next::method(@tags, @_);
}

=head2 METHODS

=head2 load

 Title   : load
 Usage   : $count = $parser->load("data.txt")
 Function: process and load data from all input files
 Returns : list of new biodatas for the sequences
 Args    : list of accessions, file names or file handles

=cut

sub load {
    my ($self) = shift;
    carp "Parsing GI list(s) and selecting accessions for download...\n" if ($self->schema->debug);
    my $term = $self->tags->{'gi file'};

    # Parse GIs
    my @input = ();
    foreach my $input (@_) {
	$input =~ s/^\s+//;
	$input =~ s/\s+.*$//;
	next unless (length $input);

	# Input is not a file or filehandle: stack it to list
	unless (UNIVERSAL::isa($input,"GLOB") || -e $input || ( ! -t STDIN && $input eq "-")) {
	    push(@input, $input);
	    next;
	}

	# Register each input filename as annotation under tag 'gi file' 
	my $fh;
	if (!UNIVERSAL::isa($input,"GLOB")) {
	    open($fh,"<$input");
	    my $file = $input eq "-" ? "standard input (pipe)" : $input;
	    foreach my $set (@{ $self->schema->datasets }) {
		my $rank = $set->search_related("dataset_attribute_values", { term_id => $term->id })->get_column("rank")->max;
		$set->find_or_create_related("dataset_attribute_values",{ term => $term, value => $file, rank  => defined $rank ? ++$rank : 0 });
	    }
	}

	# Parse input file
	while (<$fh>) {
	    chomp;
	    s/^\s+//;
	    s/\s+.*$//;
	    push(@input, $_) if (length $_);
	}
	close($fh);
    }

    # Retrieve sequences and submit to Rotifer::DBIC::AnnotationDB::Parser::fasta 
    my $nr = $self->redundant ? "F" : "T";
    my ($tmpfh, $tmpname) = ttempfile();
    id2fasta({ db => $self->blastdb, output => $tmpfh, fastacmd_opts => { "-t" => $nr, "-c" => "T" } }, @input);
    close($tmpfh);
    return $self->next::method($tmpname);
}

__PACKAGE__->meta->make_immutable;
1;
