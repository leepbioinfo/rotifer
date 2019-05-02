=head1 NAME

Rotifer::DBIC::AnnotationDB::Result::Biodata::reference - Bio::Annotatable::Reference compatible biodata

=head1 DESCRIPTION

This module emulates a Bio::Annotation::Reference modules for biodata entries.

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::DBIC::AnnotationDB::Result::Biodata::reference;

use strict;
use warnings;
use Carp::Clan qr/^Rotifer::DBIC::AnnotationDB/;
use base qw/Rotifer::DBIC::AnnotationDB::Result::Biodata/;

# This line is required!!!!
__PACKAGE__->table("biodata");

=head2 authors

 Title   : authors
 Usage   : $seq = $rs->authors()
 Function: return 1
 Returns : integer
 Args    : none

=cut

sub authors {
    return 1;
}

# Finalize Moose
1;
