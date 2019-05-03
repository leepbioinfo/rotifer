# POD documentation - main docs before the code

=head1 NAME

Rotifer::DBIC::AnnotationDB::Parser::fasta - load sequences from FASTA files 

=head1 SYNOPSIS

  # Using this parser with rannotationDB command line app
  #
  # Note: this will automatically create a dataset
  #       and add all sequences to it

  rannotationDB -a load -if fasta seqs.fasta

  # Same rannotationDB thing from a pipe

  id2fasta te | rannotationDB -a load -if fasta

  # Creating a new parser
  use Rotifer::DBIC::AnnotationDB::Parser;
  my $parser = Rotifer::DBIC::AnnotationDB::Parser->create("fasta");
  $parser->load(@ARGV);

=head1 DESCRIPTION

Rotifer::DBIC::AnnotationDB::Parser::fasta will parse and load sequences
from a list of FASTA files to a Rotifer::DBIC database. 

If the FASTA header conforms to the standards of NCBI FASTA
headers (i.e. gi|number|db|accession), it will be parsed for
NCBI identifiers and the accession numbers will be used to
track the revision history of the sequences.

By default, this parser downloads taxonomic information from
the NCBI Taxonomy Server.

=head2 Note

This class implements required methods and consumes basic methods and
attributes from AnotationDB's role for parsers. 

See Rotifer::DBIC::AnnotationDB::Role::ParserRole for detais.

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::DBIC::AnnotationDB::Parser::fasta;

use strict;
use warnings;
use autodie qw(:all);
use Bio::SeqIO;
use Carp::Clan qr/^Rotifer::DBIC/;
use DateTime;
use Moose;
use Scalar::Util qw(blessed);
use Rotifer::Utils qw(nr2ids);
use Rotifer::DB::NCBI;
use Rotifer::DB::NCBI::Taxonomy;

with ('Rotifer::DBIC::AnnotationDB::Role::BioseqPropsRole',
      'Rotifer::DBIC::AnnotationDB::Role::ParserRole' => {
	  '-excludes' => '_default_tags',
      });
    

=head2 ATTRIBUTES / ACCESSORS

This section list attributes not defined by the fundamental AnnotationDB parser's role.

=head2 download_taxonomy

 Usage   : Rotifer::DBIC::AnnotationDB::Parser->new(download_taxonomy => \%opts)
 Function: whether to download taxonomic information
 Value   : boolean
 Default : 1

=cut

has 'download_taxonomy'   => (
    is       => "rw",
    isa      => "Bool",
    default  => 1,
    );

=head2 taxopts

 Usage   : Rotifer::DBIC::AnnotationDB::Parser->new(taxopts => \%opts)
 Function: set options for download of taxonomic data
 Value   : hash reference

=cut

has 'taxopts'   => (
    is       => "rw",
    isa      => "HashRef",
    builder  => '_default_taxopts',
    required => 0
    );

sub _default_taxopts {
    return { "-update" => 0 };
}

=head2 Overriden attributes

=head2 tags

 Usage   : $parser->tags
 Function: parser annotation tags
 Returns : Rotifer::DBIC::AnnotationDB::Result::Term
 Builder : _default_tags

 The following tags will be added to the default vocabulary
 unless present in one of the other ontologies listed by
 $parser->parser_ontology (see ATTRIBUTES above):

 fasta : sequences in FASTA format

=cut

sub _default_tags {
    my $self = shift;
    return $self->_hash2tags($self->parser_ontology,
			     { name => 'fasta',
			       definition => "Sequences in FASTA format.",
			     }, @_);
}

=head2 METHODS

=head2 load

 Title   : load
 Usage   : @biodata = $parser->load("data.fa")
 Function: process and load data from all input files
 Returns : array of Rotifer::DBIC::AnnotationDB::Result::Biodata
 Args    : list of file names and/or file handlers

=cut

sub load {
    my ($self, @input) = @_;
#    my $davrs = $self->schema->resultset("DatasetAttributeValue");

    # Parsing files
    my @seqs  = ();
    my @biodata = ();
    foreach my $file (@input) {
	my $flag = UNIVERSAL::isa($file,"GLOB") ? "-fh" : "-file";
	my $parser = Bio::SeqIO->new(-format => "fasta", $flag => $file);
	while (my $seq = $parser->next_seq) {
	    push(@seqs, $seq);
	    if (scalar(@seqs) == $self->batch_size) {
		push(@biodata, $self->process_sequence_objects(@seqs));
		@seqs = ();
	    }
	}
	$parser->close;
	push(@biodata, $self->process_sequence_objects(@seqs)) if (scalar @seqs);

	# Associate input with dataset unless $file is a file handle
	if (!UNIVERSAL::isa($file,"GLOB") && scalar @biodata) {
	    $file = "standard input (pipe)" if ($file eq "-");
	    foreach my $set (@{ $self->schema->datasets }) {
		my $rank = $set->search_related("dataset_attribute_values", { term_id => $self->tags->{fasta}->id })->get_column("rank")->max;
		$set->create_related("dataset_attribute_values",
				     { term => $self->tags->{fasta}, value => $file, rank => defined $rank ? ++$rank : 0 });
	    }
	}
    }

    return @biodata;
}

=head2 process_sequence_objects

 Title   : process_sequence_objects
 Usage   : @biodata = $self->process_sequence_objects(@objs)
 Function: transfer data from Bio::Seq to Rotifer::DBIC
 Returns : integer (number of sequences loaded)
 Args    : array of Bio::SeqI

 Note    : this method is used by load()

=cut

sub process_sequence_objects {
    my ($self, @seqs) = @_;

    # Bio::Seq -> Rotifer::DBIC templates
    my @biodata = $self->prepare_biodata(@seqs);

    # Download taxonomy
    $self->schema->resultset("Biodata")
	->retrieve_all_organisms($self->taxopts, @biodata)
	if ($self->download_taxonomy);

    return @biodata;
}

=head2 prepare_biodata

 Title   : prepare_biodata
 Usage   : @biodata = $self->prepare_biodata(@objs)
 Function: transfer data from Bio::Seq to Rotifer::DBIC
 Returns : integer (number of sequences loaded)
 Args    : array of Bio::SeqI

 Note    : this method is used by process_sequence_objects()

=cut

sub prepare_biodata {
    my ($self, @seqs) = @_;
    my $xrs  = $self->schema->resultset('Dbxref');
    my $srs  = $self->schema->resultset("Biosequence");
    my $brs  = $self->schema->resultset("Biodata");
    carp "Parsing/loading ".scalar(@seqs)." FASTA sequences...\n" if ($self->schema->debug);

    # Process eqch Bio::Seq object from a fasta file: sequence, identifiers and description
    my @data   = (); # biodata stack
    my %loaded = (); # biodata registry
    foreach my $seq (@seqs) {

	# %byaccgroup contains lists of identifiers (hashes), key: accession number 
	my %byaccgroup = process_bioseq_identifiers($seq);

	# Add/retreive raw sequence text from/to database
	my $bioseq = $srs->find_or_create({ term => $self->sequence_type, residues => $seq->seq },{ cache => 1 });

	# Create/find a biodata for each accgroup
	foreach my $accgroup (keys %byaccgroup) {
	    # Just load old data unless updating taxonomy or sequence
	    my $old = $brs->find({ identifier => $accgroup, term => $self->sequence_type, rank => 0 },{ cache => 1 });
	    if (defined $old && !$self->taxopts->{"-update"}) { # If updating is not mandatory...
		if ($old->has_column_loaded("organism_id") && defined $old->get_column("organism_id")) { # If organism is set...
		    if ($old->has_column_loaded("biosequence_id") && $old->get_column("biosequence_id") eq $bioseq->id) { # If sequence is not new...
			push(@data, $old);
			next;
		    }
		}
	    }

	    # Create biodata
	    #
	    # Avoiding find_or_create() since new and update 
	    # automatically associate datasets for biodata
	    # dbxrefs, biosequences and biodata_relationship
	    my $biodata = $brs->update_or_create(
		{
		    identifier     => $accgroup,
		    term           => $self->sequence_type,
		    biosequence    => $bioseq,
		    description    => $byaccgroup{$accgroup}->[0]->{description},
		    is_rawsequence => 0,
		    is_obsolete    => 0,
		    rank           => 0,
		},{ cache => 1 });

	    # Insert dbxrefs
	    foreach my $hashref (@{$byaccgroup{$accgroup}}) {
		delete $hashref->{description}; # No such column in dbxref
		delete $hashref->{organism};    # No such column in dbxref
		my $dbxref = $xrs->update_or_create($hashref);
		$dbxref->find_or_create_related("biodata_dbxrefs"    , { biodata     => $biodata },{ cache => 1 });
		$dbxref->find_or_create_related("biosequence_dbxrefs", { biosequence => $bioseq  },{ cache => 1 });
	    }

	    # Stack biodata
	    push(@data, $biodata) unless (exists $loaded{$biodata->id});
	    $loaded{$biodata->id} = 1;
	} # foreach my $accgroup (keys %byaccgroup)
    } # foreach my $seq (@seqs)

    return @data;
}

=head2 process_bioseq_identifiers

 Title   : process_bioseq_identifiers
 Usage   : %byaccgroup = $self->process_bioseq_identifiers($seq)
 Function: parse fasta sequence identifiers
 Returns : hash of array of hashes
 Args    : Bio::SeqI complaint object

=cut

sub process_bioseq_identifiers {
    my $seq = shift;

    # nr2ids is something of a crappy parser... 
    # It returns a flat array of hashes like
    #
    # { accession => 123, acctype => "GI",        version => 0, dbname => "D", description => "A", accgroup => 123 }
    # { accession => "A", acctype => "ACCESSION", version => 0, dbname => "D", description => "A", accgroup => 123 }
    # { accession => "B", acctype => "LOCUS",     version => 0, dbname => "D", description => "A", accgroup => 123 }
    #
    # Note: when present in a group, GIs are used as accgroup

    # Parse fasta header, group hashes by accgroup
    my $header = $seq->primary_id || $seq->display_id || $seq->id || $seq->display_name;
    $header .= " ".$seq->description if ($seq->can("description") && defined $seq->description);
    my %bygroup = (); # A set of equivalent identifiers (GI, accession number and locus)
    foreach my $hashref (nr2ids($header)) { # Parse fasta header
	push(@{$bygroup{ $hashref->{accgroup} }}, $hashref);
    }

    return %bygroup;
}

__PACKAGE__->meta->make_immutable;
1;
