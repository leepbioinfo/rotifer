=head1 NAME

Rotifer::DBIC::AnnotationDB::ResultSet::Biosequence - biosequence results

=head1 DESCRIPTION
 
This modules extends DBIx::Class::ResultSet for the biosequence table.

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::DBIC::AnnotationDB::ResultSet::Biosequence;

use strict;
use warnings;
use Carp::Clan qr/^Rotifer::DBIC/;
use Digest::SHA qw(sha256_hex);
use base qw(DBIx::Class::ResultSet);

=head2 find

 Title   : find
 Usage   : $seq = $rs->find(%attrs)
 Function: find() that replaces residues by checksum 
 Returns : DBIx::Class::Row or undef
 Args    : hash reference

=cut

sub find {
    my ($self, $hash) = @_;
    $hash->{checksum} = sha256_hex($hash->{residues}) if (ref $hash && exists $hash->{residues});
    return $self->next::method($hash);
}

1;

