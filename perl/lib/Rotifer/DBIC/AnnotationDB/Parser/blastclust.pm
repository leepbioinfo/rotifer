# POD documentation - main docs before the code

=head1 NAME

Rotifer::DBIC::AnnotationDB::Parser::blastclust - load blastclust output

=head1 SYNOPSIS

  # Using this parser with rannotationDB command line app

  blastclust -S 0.3 -L 0.2 seqs.fa | rannotationDB -if blastclust

  # Creating a new parser
  use Rotifer::DBIC::AnnotationDB::Parser;
  my $parser = Rotifer::DBIC::AnnotationDB::Parser->create("blastclust");
  $parser->load(@ARGV);

=head1 DESCRIPTION

This module parsers the raw output of the blastclust program and loads
sequence clusters and, if necessary, the sequence themselves into
an instance of a Rotifer::DBIC::AnnotationDB database. 

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::DBIC::AnnotationDB::Parser::blastclust;

use strict;
use warnings;
use autodie qw(:all);
use Carp::Clan qr/^Rotifer::DBIC::AnnotationDB/;
use Moose;
extends 'Rotifer::DBIC::AnnotationDB::Parser::gi';

=head2 ATTRIBUTES

=head2 tags

 Usage   : $parser->tags
 Function: annotation tags used by this parser
 Returns : Rotifer::DBIC::AnnotationDB::Result::Term
 Builder : _default_tags

 This module adds the following tags to the parser vocabulary

  blastclust : sequence file with lots of annotations 

=cut

sub _default_tags {
    my ($self) = shift;

    my @tags = 
	(
	 {
	     name       => 'blastclust',
	     definition => "Output of NCBI's blastclust program.",
	 },
	);

    return $self->next::method(@tags, @_);
}

=head2 METHODS

=head2 load

 Title   : load
 Usage   : $full_path = $parser->load("blastclust.out")
 Function: process and load blastclust output
 Returns : 
 Args    : list of file names

=cut

sub load {
    my ($self, @input) = @_;
    my $term = $self->tags->{"genbank file"};

    # Parsing files
    my $sep     = " +";
    my %cluster = ();
    my %seqid2cluster = ();
    foreach my $input (@input) {
	my $fh;
	if (!UNIVERSAL::isa($input,"GLOB")) {
	    open($fh,"<$input");
	    my $file = $input eq "-" ? "standard input (pipe)" : $input;
	    my $term = $self->tags->{'blastclust'};
	    foreach my $set (@{ $self->schema->datasets }) {
		my $rank = $set->search_related("dataset_attribute_values", { term_id => $term->id })->get_column("rank")->max;
		$set->find_or_create_related("dataset_attribute_values",{ term => $term, value => $file, rank  => defined $rank ? ++$rank : 0 });
	    }
	}

	my $i = 0;
	open($fh, "<$input") unless UNIVERSAL::isa($fh,"GLOB");
	while (<$fh>) {
	    chomp;
	    next if /Start clustering of \d+ queries/;
	    my @F = split(/$sep/);
	    map { 
		push(@{ $cluster{$i} },$_);
		$seqid2cluster{$_} = $i;
	    } @F;
	    $i++;
	}
    }

    # Submit sequence IDs to Rotifer::DBIC::AnnotationDB::Parser::gi 
    my @biodata = $self->next::method(keys %seqid2cluster);

    # Iterate over clusters
    foreach my $id (keys %cluster) {
    }
}

__PACKAGE__->meta->make_immutable;
1;
