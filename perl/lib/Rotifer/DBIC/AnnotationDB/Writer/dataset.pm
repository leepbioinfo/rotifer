# POD documentation - main docs before the code

=head1 NAME

Rotifer::DBIC::AnnotationDB::Writer::datasets - find/print datasets 

=head1 SYNOPSIS

  # Using this parser with rannotationDB command line app
  #
  # Note: this will automatically create a dataset
  #       and add all sequences to it

  rannotationDB -of datasets -if gi 123456

  # Creating a new writer
  use Rotifer::DBIC::AnnotationDB::Writer;
  my $parser = Rotifer::DBIC::AnnotationDB::Writer->create("datasets");
  $parser->write(@ARGV);

=head1 DESCRIPTION

Search and print a table of datasets. Queries may be performed using any
valid identifier (dataset, sequence, etc.).

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

package Rotifer::DBIC::AnnotationDB::Writer::datasets;

use strict;
use warnings;
use autodie qw(:all);
use Bio::SeqIO;
use Carp::Clan qr/^Rotifer::DBIC/;
use Moose;
use Scalar::Util qw(blessed);
use Rotifer::Utils qw(ids2nr);
with qw(
 Rotifer::DBIC::AnnotationDB::Role::WriterRole
);

=head2 ATTRIBUTES / ACCESSORS

This section list attributes not defined by the fundamental AnnotationDB parser's role.

=head2 supported_input_formats

 Usage   : $writer->supported_input_formats("fasta");
 Function: list of input formats that might be used or
           from which we can extract queries to search for
 Returns : hash reference
 Value   : none, attribute is read only

 Supported input formats:

 ids
 gi
 arch

=cut

has '+supported_input_formats' => (
    default => sub {
	{
	    'Rotifer::DBIC::AnnotationDB::Result::Biodata' => 1,
	    ids  => 1,
	    gi   => 1,
	    arch => 1,
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

           Supported search methods: input, all and dataset

           Refer to rannotationDB documentation for details
           on each search method.

           Supported input formats:

           Rotifer::DBIC::AnnotationDB::Result::Biodata
           gi
           arch

=cut

sub write {
    my $self = shift;

    # Search
    my @ids = $self->parse_input(@_);
    my @all = $self->search(@ids);

    # Post-process biodata
    my %seqid = ();
    foreach my $b (@all) {
	# Add organism name to description
	if ($b->has_column_loaded("organism_id") && defined $b->get_column("organism_id")) {
	    my $orgname = $b->organism->name;
	    $b->description($b->description . " [$orgname]") if (defined $orgname);
	}

	# Merge redundant entries
	if ($b->has_column_loaded("biosequence_id") && defined $b->get_column("biosequence_id")) {
	    my $sid = $b->biosequence->id;
	    if ($self->redundant && exists $seqid{$sid}) {
		my $desc = join(" ",$seqid{$sid}->description, ">".$b->display_id, $b->description);
		$seqid{$sid}->description($desc);
		next;
	    }
	    $seqid{$sid} = $b;
	}
    }

    # Use Bio::SeqIO
    my $status = scalar(keys %seqid);
    my @opts = ("-format" => "fasta");
    unshift(@opts, "-fh", \*STDOUT) unless (exists $self->seqio_options->{"-fh"} ||
					    exists $self->seqio_options->{"-file"});
    my $io = Bio::SeqIO->new(@opts,%{ $self->seqio_options });
    foreach my $b (@seqid{ sort { $a <=> $b } keys %seqid}) {
	$io->width($b->length) if ($self->clean);
	$status = $io->write_seq($b);
    }

    return $status;
}

sub search {
    my ($self, @ids) = @_;
    my $brs = $self->schema->resultset("Datasets");
    my $batch_size = 100;
    my %bid = ();

    # Search by dataset
    my @datasets = @{ $self->schema->datasets };
    if (scalar @datasets) {
	my $start = 0;
	while ($start <= $#datasets) {
	    my $end = $start + $batch_size - 1;
	    $end = $#datasets if ($end > $#datasets);
	    my $where = {
		"me.term_id"   => $self->sequence_type->id,
		"dataset.name" => {
		    "-in" => [ map { $_->name } @datasets[$start .. $end] ],
		},
	    };
	    use Data::Dump qw(dump); carp dump($where,$attr),"\n" if ($self->schema->debug > 1);
	    map { $bid{$_->id} = $_ } $brs->search($where,$attr)->all;
	    $start = $end + 1;
	}
    } else {
	map { $bid{$_->id} = $_ } $brs->search(undef,$attr)->all;
    }

    # Search by dbxref
    $attr = { distinct => 1, join => [ { biodata_dbxrefs => "dbxref" } ]};
    if ($self->redundant) {
	$attr->{join} = [ { biosequence => { biodatas => $attr->{join}[0] } } ];
    }
    $start = 0;
    while ($start <= $#ids) {
	my $end = $start + $batch_size - 1;
	$end = $#ids if ($end > $#ids);
	my $where = {
	    "me.term_id"       => $self->sequence_type->id,
	    "dbxref.accession" => { 
		"-in" => [ @ids[$start..$end] ],
	    },
	};
	use Data::Dump qw(dump); carp dump($where,$attr),"\n" if ($self->schema->debug > 1);
	map { $bid{$_->id} = $_ if (!exists $bid{$_->id}) } $brs->search($where,$attr)->all;
	$start = $end + 1;
    }

    return sort { $b->identifier <=> $a->identifier } values %bid;
}

sub parse_input {
    my $self = shift;

    # Pre-process input
    my @data = ();
    foreach my $input (@_) {
	if (UNIVERSAL::isa($input,"GLOB") || -e $input || ( ! -t STDIN && $input eq "-")) {
	    open(my $fh,"<$input") unless (UNIVERSAL::isa($input,"GLOB"));
	    while (<$fh>) {
		chomp;
		s/^\s+//;   # Remove leading white spaces
		s/\s+.*$//; # Remove everything after the first white space
		push(@data, $_) if (length $_);
	    }
	    close($fh) unless (UNIVERSAL::isa($input,"GLOB"));
	}

	else {
	    $input =~ s/^\s+//;
	    $input =~ s/\s+$//;
	    push(@data, $input);
	}
    }

    return @data;
}

__PACKAGE__->meta->make_immutable;
1;
