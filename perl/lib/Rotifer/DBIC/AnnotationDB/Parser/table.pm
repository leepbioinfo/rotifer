# POD documentation - main docs before the code

=head1 NAME

Rotifer::DBIC::AnnotationDB::Parser::table - 
 flexible parser for text tables (planned)

=head1 SYNOPSIS

  # Using this parser through rannotationDB
  rannotationDB -if table table.tsv

  # A slightly more sophisticated example
  rannotationDB -f table --input_args delimiter=" : " table.tsv

  # If you would like to program:
  #
  # Creating a new parser
  use Rotifer::DBIC::AnnotationDB::Parser;
  my $parser = Rotifer::DBIC::AnnotationDB::Parser->parser(format => "table");
  $parser->load(%options, @ARGV);

=head1 DESCRIPTION

Rotifer::DBIC::AnnotationDB::Parser::table parses text tables and 
inserts or updates the parsed data. It requires the right column
names to correctly map data to a Rotifer::DBIC::AnnotationDB
instance.

=head2 Column headers in your input table 

How do you name your columns in the input table to correctly store your data?

Column names should be formatted as 'table name.column name', 
where 'table' is the name of a table in Rotifer::DBIC::AnnotationDB
and 'column' is the name of a column in that table.

A detailed description of all tables and columns in AnnotationDB can
be found in rotifer's share/doc directory.

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 DEPENDENCIES

=over

=item Moose

=item Rotifer::DBIC

=back

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::DBIC::AnnotationDB::Parser::table;

use strict;
use warnings;
use namespace::autoclean;
use autodie qw(:all);
use Moose;

=head1 METHODS

=head2 load

 Title   : load
 Usage   : $full_path = $parser->load("data.txt")
 Function: process and load data from all input files
 Returns : 
 Args    : list of file names

=cut

sub load {
    my ($self, @input) = @_;

    # Loading all tables
    foreach my $file (@input) {
#	$self->_parse_table($file);
    }

    return 1;
}

sub _parse_table {
    my ($self, $file) = @_;

    my $row_number = 0;
    my $sep = $self->delimiter;
    my $processed_header = undef;
    open(my $th,"<$file");
    while (<$th>) {
	chomp;
	$row_number++; # Set the current row number
	my @row = split(/$sep/, $_, $self->number_of_columns);
	@row = map { s/^\s+//; s/\s+$//; $_ } @row if $self->unpad;

	# Load header from table unless set by the user
	if (!defined $processed_header) {
	    $processed_header = $self->_parse_header(scalar @{$self->header} ? @{$self->header} : @row);
	    next;
	}

	# Check consistency
	if ($#row > $#{$self->header}) {
	    warn "WARNING: row $row_number exceeds the number of columns named in the header of table $file. Unnamed columns will NOT be loaded!!";
	}

    } # while (<$th>)
    close($th);
} # sub _load_table

# Process each file's header and plan queries
# 
sub _parse_header {
    my ($self, @header) = @_;

    my $processed = {};
    for (my $i=0; $i<=$#header; $i++) {
	my ($table, $column) = split(/\./,$header[$i]);
	$processed->{$table}{$i} = $column;
    }

    return $processed;
}

__PACKAGE__->meta->make_immutable;
1;
