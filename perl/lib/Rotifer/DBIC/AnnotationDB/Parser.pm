# POD documentation - main docs before the code

=head1 NAME

Rotifer::DBIC::AnnotationDB::Parser - Rotifer::DBIC::AnnotationDB parser factory

=head1 SYNOPSIS

  # Creating a new parser
  use Rotifer::DBIC::AnnotationDB::Parser;
  my $parser = Rotifer::DBIC::AnnotationDB::Parser->create("fasta", \%opts);
  $parser->load(@ARGV);

=head1 DESCRIPTION

Rotifer::DBIC::AnnotationDB::Parser chooses and instantiates parsers for
different file formats. These parsers load, insert and
update data in a Rotifer::DBIC database.

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 DEPENDENCIES

=over

=item MooseX::AbstractFactory

=back

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::DBIC::AnnotationDB::Parser;

use strict;
use warnings;
use MooseX::AbstractFactory;
use Rotifer::Utils qw(describe_subclasses);

# Set my subclass namespace
implementation_class_via sub { 
    my $parser = shift;
    if ($parser eq "help") {
	print describe_subclasses(__PACKAGE__);
	exit 1;
    }
    return __PACKAGE__."::$parser";
};

__PACKAGE__->meta->make_immutable;
1;
