# POD documentation - main docs before the code

=head1 NAME

Rotifer::DBIC::AnnotationDB::Search - Rotifer::DBIC::AnnotationDB search factory

=head1 SYNOPSIS

  # Creating a new search
  use Rotifer::DBIC::AnnotationDB::Search;
  my $search = Rotifer::DBIC::AnnotationDB::Search->create("fasta", \%opts);
  $search->load(@ARGV);

=head1 DESCRIPTION

Rotifer::DBIC::AnnotationDB::Search chooses and instantiates 
pre-defined searches for a Rotifer::DBIC database. 

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

package Rotifer::DBIC::AnnotationDB::Search;

use strict;
use warnings;
use MooseX::AbstractFactory;
use Rotifer::Utils qw(describe_subclasses);

# Set my subclass namespace
implementation_class_via sub { 
    my $search = shift;
    if ($search eq "help") {
	print describe_subclasses(__PACKAGE__);
	exit 1;
    }
    return __PACKAGE__."::$search";
};

__PACKAGE__->meta->make_immutable;
1;
