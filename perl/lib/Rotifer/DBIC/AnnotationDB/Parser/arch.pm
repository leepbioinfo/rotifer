# POD documentation - main docs before the code

=head1 NAME

Rotifer::DBIC::AnnotationDB::Parser::arch - load domain architectures

=head1 SYNOPSIS

  # Using this parser with rannotationDB command line app
  #
  # Note: this will automatically create a biograph
  #       and add all sequences to it

  rannotationDB -a load -if arch te.arch

  # Same thing from a pipe

  domain2architecture hits.tsv | rannotationDB -a load -if arch

  # Creating a new parser
  use Rotifer::DBIC::AnnotationDB::Parser;
  my $parser = Rotifer::DBIC::AnnotationDB::Parser->create("arch");
  $parser->load(@ARGV);

=head1 DESCRIPTION

Rotifer::DBIC::AnnotationDB::Parser::arch parses and load architecture tables generated
by domain2architecture and compatible programs.

Architetures must be in a format equivalent to the output of domain2architecture.
The minimum requirements are two columns:

=over

=item (1) The first column must contain the sequence identifiers (usually GI or fake GIs)

=item (2) The second column must contain architectures as seen in domain2architecture's output

=back

Additional columns produced by domain2architecture are supported and include

=over

=item (3) A third column with coordinates along the sequences identified in the first column

=item (4) The e-value scores for each match (loaded into column score, table biodata_relationship)

=item (5) A fifth column is also supported and may contain coordinates over the length of
profiles or queries used to identify each domain.

=back

=head1 EXTENDS

 Rotifer::DBIC::AnnotationDB::Parser::gi

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::DBIC::AnnotationDB::Parser::arch;

use strict;
use warnings;
use autodie qw(:all);
use Carp::Clan qr/^Rotifer::DBIC::AnnotationDB/;
use Data::UUID;
use List::MoreUtils qw(uniq);
use Moose;
use Scalar::Util qw(blessed);
use Rotifer::DBIC::AnnotationDB;
extends 'Rotifer::DBIC::AnnotationDB::Parser::gi';

=head2 ATTRIBUTES

=head2 protein_domain_ontology

 Usage   : Rotifer::DBIC::AnnotationDB::Parser->new(protein_domain_ontology => "Protein Universe")
 Function: input column separator
 Value   : string
 Builder : _default_protein_domain_ontology
 Trigger : _protein_domain_ontology_trigger

=cut

has 'protein_domain_ontology' => (
    is       => 'rw',
    isa      => 'Str|HashRef|Object',
    lazy     => 1,
    builder  => '_default_protein_domain_ontology',
    trigger  => \&_protein_domain_ontology_trigger,
    );

sub _default_protein_domain_ontology {
    return shift->schema->resultset('Ontology')->
	find_or_create(name => "protein domain ontology (sensu $ENV{USER})");
}

sub _protein_domain_ontology_trigger {
    my ($self, $new, $old) = @_;
    unless (blessed $new && $new->isa("Rotifer::DBIC::AnnotationDB::Result::Ontology")) {
	my $hash = { name => $new } unless UNIVERSAL::isa($new,"HASH");
	my $onto = $self->schema->resultset('Ontology')->find_or_create($hash);
	$self->protein_domain_ontology($onto);
    }
}

=head2 protein_domain_tags

 Usage   : Rotifer::DBIC::AnnotationDB::Parser->new(protein_domain_tags => \@objs)
 Function: find/create base terms for protein domain annotation
 Value   : array reference of Terms
 Builder : _default_protein_domain_tags

=cut

has 'protein_domain_tags' => (
    is       => "ro", 
    isa      => "HashRef[Object]",
    lazy     => 1,
    builder  => '_default_protein_domain_tags',
    init_arg => undef,
    );

sub _default_protein_domain_tags {
    my $self = shift;

    my @tags = 
	(
	 {
	     name       => 'PART_OF',
	     definition => 'A is part of B if for every A there is some subcomponent of A named B',
	     synonyms   => [ 'part_of' ],
	 },
	 {
	     name       => "HAS_PART",
	     definition => 'Inverse of part_of',
	     synonyms   => [ 'has_part' ],
	 },
	 {
	     name       => 'protein domain architecture',
	     definition => 'A set of protein domains present in a given protein sequence.',
	 },
	 {
	     name       => 'protein domain model',
	     definition => 'A mathematical model of a collection of related protein sequences. E.g. one or more sequence(s), a PSSM or a HMM.',
	 },
	);

    return $self->_hash2tags($self->protein_domain_ontology, @tags);
}

=head2 add_domains_to_ontology

 Usage   : Parser->new(add_domains_to_ontology => 1)
 Function: if true, add domain names to the ontology (create terms)
           if false, abort and complain if domain names are not found
 Value   : boolean

=cut

has 'add_domains_to_ontology' => (
    is       => 'rw',
    isa      => 'Bool',
    default  => 1,
    );

=head2 input_delimiter

 Usage   : Rotifer::DBIC::AnnotationDB::Parser->new(input_delimiter => "\t")
 Function: input column separator
 Value   : string

=cut

has input_delimiter => (
    is      => 'rw',
    isa     => 'Str',
    default => "\t",
    );

=head2 tags

 Usage   : $parser->tags
 Function: annotation tags used by this parser
 Returns : Rotifer::DBIC::AnnotationDB::Result::Term
 Builder : _default_tags

=cut

sub _default_tags {
    my ($self) = shift;

    my @tags = 
	(
	 {
	     name       => 'architecture table',
	     definition => "Description of protein domain architectures generated by domain2architecture.",
	 },
	);

    return $self->next::method(@tags, @_);
}

=head2 METHODS

=head2 load

 Title   : load
 Usage   : $full_path = $parser->load("data.txt")
 Function: process and load data from all input files
 Returns : 
 Args    : list of file names

=cut

sub load {
    my ($self, @input) = @_;

    # Parse input
    my @archs = $self->parse_architectures(@input);

    # Process archs in batches
    my $start = 0;
    my @ret = ();
    while ($start <= $#archs) {
	# Prepare current batch
	my $end = $start + $self->batch_size - 1;
	$end = $#archs if ($end > $#archs);
	my @batch = @archs[${start} .. $end];

	# Add sequences
	my @ids = uniq map { $_->[0]{object}{identifier} } grep { $_->[0]{predicate}->name eq 'PART_OF' } @batch;
	$self->next::method(@ids); # use load() from Rotifer::DBIC::AnnotationDB::Parser::gi

	# Add domains
	my @biodata = $self->process_architectures(@batch);

	# Prepare next batch interval
	push(@ret, @biodata);
	$start = $end + 1;
    }

    return @ret;
}

=head2 parse_architectures

 Title   : parse_architectures
 Usage   : $full_path = $self->parse_architectures("data.txt")
 Function: parse architecture table to one array of arrays of hashes
 Returns : hash of arrays (input table's first column is the key) 
 Args    : list of file names

 NOTE: input must be identical to the output of domain2architecture 

=cut

sub parse_architectures {
    my $self = shift;
    carp "Parsing ".scalar(@_)." files for domain architectures..." if ($self->schema->debug);
    my $sep  = $self->input_delimiter;
    my $trs  = $self->schema->resultset('Term');
    my $brs  = $self->schema->resultset('Biodata');

    my @domains = ();
    foreach my $file (@_) {
	my @unknown = ();
	my $archFH  = $file;
	open($archFH,"<$archFH") unless UNIVERSAL::isa($archFH,"GLOB");
	while (<$archFH>) {
	    chomp;
	    next if /^\s*$/;

	    # Parse row: first column = source ID
	    my @row = split(/$sep/);
	    next unless ($#row > 0);

	    # Object to represent the architecture
	    my $object = {
		identifier     => $row[0],
		term           => $self->protein_domain_tags->{'protein domain architecture'},
		rank           => 0,
		is_obsolete    => 0,
		is_rawsequence => 1,
	    };

	    # Third column: domain coordinates
	    my @start = (); my @end = ();
	    if (defined $row[2]) {
		my @regions = split(/\s*\,\s*/,$row[2]);
		for (my $i=0; $i<=$#regions; $i++) {
		    my ($start, $end, $name) = split(/\s*[\.\&]+\s*/,$regions[$i]);
		    $start[$i] = $start;
		    $end[$i]   = $end;
		}
	    }

	    # Fourth column: e-value scores
	    my @score = ();
	    if (defined $row[3]) {
		my @score = split(/\s*\,\s*/,$row[2]);
		for (my $i=0; $i<=$#score; $i++) {
		    $score[$i] =~ s/^e-/1e-/;
		}
	    }

	    # Second column: compact description of architectures is 
	    # expected to be in the second column of the input table 
	    my $rank  = 0;
	    my $match = [];
	    my @name = grep { defined $_ } split(/[\+\[\]]+/,$row[1]);
	    for (my $i=0; $i<=$#name; $i++) {
		# Check domain names
		my $cond = { name        => $name[$i],
			     ontology_id => $self->protein_domain_ontology->id,
			     is_obsolete => 0 };
		my $term = $trs->find($cond, { cache => 1 });
		if (!defined $term) {
		    if ($self->add_domains_to_ontology) {
			$term = $trs->create($cond, { cache => 1 });
		    } else {
			push(@unknown, $name[$i]);
		    }
		}
		next if (scalar @unknown);

		# Create relationship
		my $forward = {
		    subject => {
			term => $term,
			rank => 0,
			is_rawsequence => 0,
			is_obsolete => 0,
		    },
		    predicate    => $self->protein_domain_tags->{'PART_OF'},
		    object       => $object,
		    rank         => $i,
		    start_coord  => $start[$i],
		    end_coord    => $end[$i],
		    strand       => 1,
		    score        => $score[$i],
		    is_location  => 1,
		};
		push(@$match, $forward);
	    } # for (my $i=0; $i<=$#name; $i++)
	    next unless (scalar @$match);

	    # Fifth column: planned
	    #
	    # Architecture files store the most general description of
	    # a domain name in their second column. The best model name
	    # is usualy lost but an extension of the format could allow
	    # to store the model name and coordinates on the model in
	    # the fifth column

	    if (defined $row[4]) {
		my @regions = split(/\s*\,\s*/,$row[4]);
		for (my $i=0; $i<=$#regions; $i++) {
		    my ($start, $end, $name) = split(/\s*[\.\&]+\s*/,$regions[$i]);
		    my $model = $brs->find_or_create({ identifier     => $name,
						       term           => $self->protein_domain_tags->{'protein domain model'},
						       rank           => 0,
						       is_obsolete    => 0,
						       is_rawsequence => 1,
						     });
		    my $reverse = {
			subject      => $object,
			predicate    => $self->protein_domain_tags->{'HAS_PART'},
			object       => $model,
			rank         => $i,
			start_coord  => $start,
			end_coord    => $end,
			strand       => 1,
			score        => $score[$i],
			is_location  => 1,
		    };

		    push(@$match, $reverse);
		} # for (my $i=0; $i<=$#regions; $i++)
	    } # if (defined $row[4])

	    # Add relationships
	    push(@domains, $match) if (scalar @$match);
	} # while (<$archFH>)
	close($archFH);

	die join("\n","The following domain names were not found in your ontology:",@unknown,
		 "You may add these terms yourself or enable autocreation with 'add_domains_to_ontology=1' input parser option")
	    if (scalar @unknown);

	# Register data origin
	if (!UNIVERSAL::isa($file,"GLOB")) {
	    $file = "standard input (pipe)" if ($file eq "-");
	    foreach my $dataset (@{ $self->schema->datasets }) {
		my $rank = $dataset->search_related("dataset_attribute_values",{ term_id => $self->tags->{"architecture table"}->id })->get_column("rank")->max;
		$dataset->create_related("dataset_attribute_values",{ term => $self->tags->{"architecture table"}, value => $file, rank => defined $rank ? ++$rank : 0 });
	    }
	}
    } # foreach my $file (@_)

    return @domains;
}

=head2 process_architecturess

 Title   : process_architecturess
 Usage   : $full_path = $self->process_architecturess("data.txt")
 Function: parse architecture table into an array of hashes 
 Returns : hash of arrays (input table's first column is the key) 
 Args    : list of file names

=cut

sub process_architectures {
    my ($self, @archs) = @_;
    carp "Processing ".scalar(@archs)." architectures..." if ($self->schema->debug);

    # Load ResultSets
    my $brs  = $self->schema->resultset('Biodata');
    my $brrs = $self->schema->resultset('BiodataRelationship');
    my $xrs  = $self->schema->resultset('Dbxref');

    # Process each accession
    my @data = ();
    foreach my $arch (@archs) {
	# Since architecture files have no sequences, search biosequence by identifier
	my @biodata = $xrs->get_seqs_by_id($arch->[0]{object}{identifier});
	next unless (scalar @biodata);
	die "Architecture files should use unique identifiers for sequences but " . 
	    $arch->[0]{object}{identifier} . " is associated with ".scalar(@biodata) .
	    " sequences" if (scalar @biodata > 1);
	my $bioseq = $biodata[0]->biosequence;
	next unless (defined $bioseq);
	my $dbxref = $xrs->search({ accession => $arch->[0]{object}{identifier} });

	# Update new architecture
	my $newarch = $arch->[0]{object};
	$newarch->{identifier}  = "arch:" . $bioseq->id;
	$newarch->{biosequence} = $bioseq;
	map { $_->{subject}{biosequence} = $bioseq } @$arch;

	# Retrieve old architecture
	my @oldarch = $bioseq->search_related
	    ("biodatas", {
		term_id        => $self->protein_domain_tags->{"protein domain architecture"}->id,
		is_obsolete    => 0,
		is_rawsequence => 1,
	     }, { order_by => [qw/rank/] })->all;

	# Compare old and new architecture
	if (scalar @oldarch) {
	    $newarch->{rank} = $oldarch[$#oldarch]->rank + 1;
	    if ($oldarch[$#oldarch]->is_equal_to($arch)) {
		push(@data, $oldarch[$#oldarch]);
		next;
	    }
	    map { $_->is_obsolete(1); $_->update } @oldarch; # Fix any older non-obsolete architectures
	}

	# Insert new architecture
	map {
	    $_->{subject}{identifier} = join(":", "domain", $bioseq->id, $newarch->{rank}, $_->{rank});
	    $_->{subject} = $brs->find_or_create($_->{subject}, { cache => 1 })
	} @$arch;
	$newarch = $brs->create($newarch);
	map { $brrs->create($_) } @$arch;
	map { $_->find_or_create_related("biodata_dbxrefs", { biodata => $newarch }) } $dbxref->all if (defined $dbxref);
	push(@data, $newarch);
    } # foreach my $arch (@archs)

    return @data;
}

__PACKAGE__->meta->make_immutable;
1;
