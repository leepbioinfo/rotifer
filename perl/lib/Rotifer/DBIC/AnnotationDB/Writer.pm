# POD documentation - main docs before the code

=head1 NAME

Rotifer::DBIC::AnnotationDB::Writer - Rotifer::DBIC::AnnotationDB writer factory

=head1 SYNOPSIS

  # Creating a new writer
  use Rotifer::DBIC::AnnotationDB::Writer;
  my $writer = Rotifer::DBIC::AnnotationDB::Writer->create("fasta", \%opts);
  $writer->load(@ARGV);

=head1 DESCRIPTION

Rotifer::DBIC::AnnotationDB::Writer chooses and instantiates writers for
different file formats. A writer search  a Rotifer::DBIC database based
on user defined parameters and formats the results in a suitable output
format.

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

package Rotifer::DBIC::AnnotationDB::Writer;

use strict;
use warnings;
use MooseX::AbstractFactory;
use Rotifer::Utils qw(describe_subclasses);

# Set my subclass namespace
implementation_class_via sub { 
    my $writer = shift;
    if ($writer eq "help") {
	print describe_subclasses(__PACKAGE__);
	exit 1;
    }
    return __PACKAGE__."::$writer";
};

__PACKAGE__->meta->make_immutable;
1;
