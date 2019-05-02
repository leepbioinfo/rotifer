# POD documentation - main docs before the code

=head1 NAME

Rotifer::DBIC::AnnotationDB::Writer::conflicts - 
 find proteins annotated with conflicting architectures

=head1 SYNOPSIS

  # Using this parser with rannotationDB command line app
  #
  # Note: this will automatically create a dataset
  #       and add all architectures to it

  rannotationDB -of conflicts -if arch all.arch

  # Creating a new writer
  use Rotifer::DBIC::AnnotationDB::Writer;
  my $parser = Rotifer::DBIC::AnnotationDB::Writer->create("conflicts");
  $parser->write(@ARGV);

=head1 DESCRIPTION

This writer will search a Rotifer::DBIC::AnnotationDB database for 
conflicting protein architectures of all proteins in a dataset.

Users may choose to search for conflicts using one or more ontologies
to annotate sequence domains (see attribute protein_domain_ontology). 

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

package Rotifer::DBIC::AnnotationDB::Writer::conflicts;

use strict;
use warnings;
use autodie qw(:all);
use Carp::Clan qr/^Rotifer::DBIC/;
use Moose;
use Scalar::Util qw(blessed);
use Rotifer::Utils qw(aoa2txt aoa2tsv);
with 'Rotifer::DBIC::AnnotationDB::Role::WriterRole';

=head2 ATTRIBUTES / ACCESSORS

This section list attributes not defined by the fundamental AnnotationDB parser's role.

=head2 supported_input_formats

 Usage   : $writer->supported_input_formats("fasta");
 Function: list of input formats that might be used or
           from which we can extract queries to search for
 Returns : hash reference
 Value   : none, attribute is read only

Supported formats are: arch, gi, ids, domains

=cut

has '+supported_input_formats' => (
    default => sub {
	{
	    'Rotifer::DBIC::AnnotationDB::Result::Biodata' => 1,
	    arch    => 1,
	    domains => 1,
	    gi      => 1,
	    ids     => 1,
	}
    });

=head2 METHODS

=head2 write

 Title   : write
 Usage   : @biodata = $parser->write(@args)
 Function: process input, search and write data
 Returns : true on success, false otherwise
 Args    : array of file names and/or sequence identifiers 

 Notes   : input files could be any table with sequence
           identifiers in the first column, like GI lists
           or architecture tables. Name the right type of
           ID using --input_format (see suuported format
           above)

=cut

sub write {
    my $self = shift;
    my $sql = Rotifer->sysdir('lib') . "/sql/AnnotationDB/conflicts." . lc($self->schema->storage->sqlt_type);
    open(my $fh, "<$sql");
    $sql = join("",<$fh>);
    close($fh);
    $self->schema->storage->dbh_do(
	sub {
	    my ($storage, $dbh, @cols) = @_;
	    my $ref = $dbh->selectall_arrayref($sql);
	    my %group = ();
	    foreach my $row (@$ref) {
		push(@{ $group{$row->[1]} },$row);
	    }
	    foreach my $gid (sort { scalar(@{$group{$b}}) <=> scalar(@{$group{$a}}) } keys %group) {
		my $arch = undef; my $nofDom = 0;
		foreach my $row (@{$group{$gid}}) {
		    my @d = split(/\+/, $row->[2]);
		    if (!defined $arch || $nofDom < scalar(@d) || length($arch) < length($row->[2])) {
			$arch = $row->[2];
			$nofDom = scalar(@d);
		    }
		    splice(@$row,1,1);
		}
		print "# ",scalar(@{$group{$gid}}),"; ",$arch,"\n";
		print aoa2tsv(@{$group{$gid}});
	    }
	});
}

__PACKAGE__->meta->make_immutable;
1;
