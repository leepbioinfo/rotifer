# POD documentation - main docs before the code

=head1 NAME

Rotifer::DBIC::AnnotationDB::Writer::phyloprofile -
    distributions of sequence features among taxa

=head1 SYNOPSIS

  # Default: get the number of proteins per organism

  rannotationDB -of phyloprofile -project myproj1

  # Include also the distribution of domains
  # See rotifer -doc Rotifer::DBIC::AnnotationDB::Parser::arch

  rannotationDB -of phyloprofile -project myproj1 -oa domains=1

  # Include genes (ideally, use a sequence ontology)

  rannotationDB -of phyloprofile -project myproj1 -oa genes=1

  # Creating a new writer
  use Rotifer::DBIC::AnnotationDB::Writer;
  my $parser = Rotifer::DBIC::AnnotationDB::Writer->create("phyloprofile");
  $parser->write(@ARGV);

=head1 DESCRIPTION

This writer will search a Rotifer::DBIC::AnnotationDB database and
compile count the number of sequence features for each taxa. 

=head2 Note

This class implements required methods and consumes basic methods and
attributes from AnotationDB's role for writers. 

See Rotifer::DBIC::AnnotationDB::Role::WriterRole for detais.

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::DBIC::AnnotationDB::Writer::phyloprofile;

use strict;
use warnings;
use autodie qw(:all);
use Carp::Clan qr/^Rotifer::DBIC/;
use Moose;
use Scalar::Util qw(blessed);
with 'Rotifer::DBIC::AnnotationDB::Role::WriterRole';

=head2 ATTRIBUTES / ACCESSORS

This section list attributes not defined by the fundamental AnnotationDB parser's role.

=head2 supported_input_formats

 Usage   : $writer->supported_input_formats("arch");
 Function: list of input formats that might be used or
           from which we can extract queries to search for
 Returns : hash reference
 Value   : none, attribute is read only

 Supported input formats are

  domain : terms or synonyms from the protein domain ontology
  ids    : any Dbxref accession that can be used to fetch sequences
  gi     : NCBI's GI or accession numbers
  arch   : a protein domain architecture file

=cut

has '+supported_input_formats' => (
    default => sub {
	{
	    domain => 1,
	    ids    => 1,
	    gi     => 1,
	    arch   => 1,
	}
    });

=head2 METHODS

=head2 write

 Title   : write
 Usage   : @biodata = $parser->write(@args)
 Function: process input, search and write data
 Returns : true on success, false otherwise
 Args    : array of sequence identifiers, files or biodata

 Notes   : input files could be any table with sequence
           identifiers in the first column, like GI lists
           or architecture tables.

           Supported input formats:

           domain
           ids
           gi
           arch

=cut

sub write {
    my $self = shift;
    my @ids = $self->parse_ids(@_);
    my $brs = $self->schema->resultset('Biodata');

    # Search
    my %arch = (); my %gi = ();
    foreach my $cursor ($self->search_architectures(@ids)) {
	while (my @data = $cursor->next) {
	    next if (exists $gi{$data[1]});
	    my $arch = "";
	    if (exists $arch{$data[0]}) {
		$arch = $arch{$data[0]};
	    } else {
		$arch = $brs->find($data[0])->architecture;
		$arch{$data[0]} = $arch;
	    }
	    print join("\t",$data[1],$arch),"\n";
	    $gi{$data[1]} = 1;
	}
    }

    return 1;
}

sub search_architectures {
    my ($self,@ids) = @_;
    my $brs = $self->schema->resultset("Biodata");

    # Where
    my $where = {
	"me.is_obsolete" => 0,
	"term.name"      => "protein domain architecture",
    };

    # Join
    my $join = 
	[
	 "term", 
	 { biosequence => { "biosequence_dbxrefs" => "dbxref" } }
	];

    # Restrict by datasets
    if (scalar @{ $self->schema->datasets }) {
	$join = [ @$join, { biodata_datasets => "dataset" } ];
	$where->{"dataset.name"} = { 
	    "-in" => [ map { $_->name } @{ $self->schema->datasets } ]
	}
    }

    # Restrict by domain name
    if (defined $self->input_format && $self->input_format eq "domain") {
	$join = [ @$join, { biodata_relationship_objects => { subject => "term" } } ];
    }

    # Process ids and search
    my $batch = 100;
    my $start = scalar(@ids) ? 0 : -1;
    my @cursor = ();
    while ($start <= $#ids) {
	my $cond = { %$where };

	# Restrict by ID
	if (scalar @ids) {
	    my $end = $start + $batch - 1;
	    $end = $#ids if ($end > $#ids);
	    my $col = "dbxref.accession";
	    $col = "term_2.name" if ($self->input_format eq "domain");
	    $cond->{$col} = { "-in" => [ @ids[$start..$end] ] };
	    if (exists $cond->{"dataset.name"}) {
		$cond->{"-or"} = {
		    "dataset.name" => delete $cond->{"dataset.name"},
		    $col => delete $cond->{$col},
		};
	    }
	}

	# Search
	my $attr = {
	    join     => $join,
	    columns  => [qw/me.biodata_id dbxref.accgroup/],
	    distinct => 1,
	};
	use Data::Dump qw(dump); print dump(($cond, $attr)),"\n" if ($self->schema->debug > 1);
	my $arch_rs = $brs->search($cond,$attr);

	push(@cursor, $arch_rs->cursor);
	$start += $batch;
    }

    return @cursor;
}

sub parse_ids {
    my $self = shift;

    # Pre-process input
    my %data = ();
    foreach my $input (@_) {

	# Input is a file name or file handler
	if (UNIVERSAL::isa($input,"GLOB") || -e $input || ( ! -t STDIN && $input eq "-")) {
	    open(my $fh,"<$input") unless (UNIVERSAL::isa($input,"GLOB"));
	    while (<$fh>) {
		chomp;
		s/^\s+//;   # Remove leading white spaces
		s/\s+.*$//; # Remove everything after the first white space
		next unless (length $_);
		$data{$_}++;
	    }
	    close($fh) unless (UNIVERSAL::isa($input,"GLOB"));
	}

	# Input is a string
	else {
	    $input =~ s/^\s+//;
	    $input =~ s/\s+$//;
	    $data{$input}++;
	} 

    } # foreach my $input (@_)

    # Remove repetead queries
    return keys %data;
}

__PACKAGE__->meta->make_immutable;
1;
