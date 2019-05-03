# POD documentation - main docs before the code

=head1 NAME

Rotifer::DBIC::AnnotationDB::Search::fasta - write sequences to FASTA format 

=head1 SYNOPSIS

  # Using this parser with rannotationDB command line app
  #
  # Note: this will automatically create a dataset
  #       and add all sequences to it

  rannotationDB -of fasta 123456

  # Creating a new search
  use Rotifer::DBIC::AnnotationDB::Search;
  my $parser = Rotifer::DBIC::AnnotationDB::Search->create("fasta");
  $parser->write(@ARGV);

=head1 DESCRIPTION

This search will search a Rotifer::DBIC::AnnotationDB database for biodata
matching a list of identifiers (accession numbers, GIs, etc.) and print
any matching biodata in FASTA format.

=head2 Note

This class implements required methods and consumes basic methods and
attributes from AnotationDB's role for searchs. 

See Rotifer::DBIC::AnnotationDB::Role::SearchRole for detais.

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::DBIC::AnnotationDB::Search::fasta;

use strict;
use warnings;
use autodie qw(:all);
use Bio::SeqIO;
use Carp::Clan qr/^Rotifer::DBIC/;
use Moose;
use Scalar::Util qw(blessed);
use Rotifer::Utils qw(ids2nr);
with qw(
 Rotifer::DBIC::AnnotationDB::Role::SearchRole
 Rotifer::DBIC::AnnotationDB::Role::BioseqPropsRole
);

=head2 ATTRIBUTES / ACCESSORS

This section list attributes not defined by the fundamental AnnotationDB parser's role.

=head2 supported_input_formats

 Usage   : $search->supported_input_formats("fasta");
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

=head2 redundant

 Usage   : Rotifer::DBIC::AnnotationDB::Search->new(redundant => 1);
 Function: return all identical sequences matching your query
 Returns : boolean
 Value   : boolean
 Default : false (0)

=cut

has 'redundant' => (
    is      => "rw", 
    isa     => "Bool",
    default => 0,
    );

=head2 seqio_options

 Usage   : Rotifer::DBIC::AnnotationDB::Search->new(seqio_options => { "-file" => "a.fa" })
 Function: get/set Bio::SeqIO options
 Returns : hash reference
 Value   : hash
 Builder : _default_seqio_options
 Trigger : _seqio_options_trigger

=cut

has 'seqio_options' => (
    is      => "rw", 
    isa     => "HashRef",
    lazy    => 1,
    builder => '_default_seqio_options',
    trigger => \&_seqio_options_trigger,
    );

sub _default_seqio_options {
    return {};
}

sub _seqio_options_trigger {
    my ($self, $new, $old) = @_;
    foreach my $key (keys %$new) {
	next if ($key =~ /^-+/);
	my $value = delete $new->{"$key"};
	$new->{"-$key"} = $value;
    }
}

=head2 METHODS

=head2 search

 Title   : search
 Usage   : @biodata = $parser->search(@args)
 Function: process input and search data
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

sub search {
    my ($self, @ids) = @_;
    my $brs = $self->schema->resultset("Biodata");
    my $batch_size = 100;
    my %bid = ();

    # Search by dataset
    my $attr = { distinct => 1, join => [ { biodata_datasets => "dataset" } ]};
    if ($self->redundant) {
	$attr->{join} = [ { biosequence => { biodatas => $attr->{join}[0] } } ];
    }
    my $start = 0;
    my @datasets = @{ $self->schema->datasets };
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
