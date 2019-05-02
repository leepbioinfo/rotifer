=head1 NAME

Rotifer::DBIC::AnnotationDB::ResultSet::Biodata - biodata resultset

=head1 DESCRIPTION
 
This modules extends DBIx::Class::ResultSet for the table biodata.

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::DBIC::AnnotationDB::ResultSet::Dbxref;

use strict;
use warnings;
use Carp::Clan qr/^Rotifer::DBIC/;
use Moose;
use MooseX::NonMoose;
use MooseX::MarkAsMethods autoclean => 1;
extends 'DBIx::Class::ResultSet';

=head1 Utilities

=head2 get_seqs_by_id

 Title   : get_seqs_by_id
 Usage   : $biodata_rs = $self->get_seqs_by_id($identifier)
 Function: retrieve biodata
 Returns : array of Rotifer::DBIC::AnnotationDB::Result::Biodata
 Args    : string

=cut

sub get_seqs_by_id {
    my $self = shift;
    my $batch_size = 100;
    my $start = 0;
    my $end   = $batch_size - 1;
    my @data = ();
    while ($start <= $#_) {
	$end = $#_ if ($end > $#_);
	my @batch = @_[$start..$end];
	push(@data, $self
	     ->search_related("biodata_dbxrefs")
	     ->search_related("biodata",
			      {
				  accession => { '-in' => \@batch },
				  '-or' => {
				      "term.name"             => { 'like' => [qw( dna protein %rna )] },
				      "term_synonyms.synonym" => { 'like' => [qw( dna protein %rna )] },
				  },
			      },{
				  join     => { term => "term_synonyms" },
				  distinct => 1,
			      })
	     ->all);
	$start = $end + 1;
	$end  += $batch_size;
    }
    return @data;
}

__PACKAGE__->meta->make_immutable();
1;
