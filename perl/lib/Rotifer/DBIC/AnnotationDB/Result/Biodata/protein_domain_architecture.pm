=head1 NAME

Rotifer::DBIC::AnnotationDB::Result::Biodata::protein_domain_architecture - 
 biodata extension for protein domains

=head1 DESCRIPTION

This module implemets helper methods for Biodata objects to deal with
the specifics of protein domain architectures.

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::DBIC::AnnotationDB::Result::Biodata::protein_domain_architecture;

use strict;
use warnings;
use Carp::Clan qr/^Rotifer::DBIC::AnnotationDB/;
use Scalar::Util qw(blessed);
use base qw/Rotifer::DBIC::AnnotationDB::Result::Biodata/;

# This line is required!!!!
__PACKAGE__->table("biodata");

=head1 ATTRIBUTES

=head1 METHODS

=head2 architecture

 Title   : architecture
 Usage   : $str = $biodata->architecture()
 Function: compact representation of a protein domain architecture
 Returns : string

=cut

sub architecture {
    my ($self) = @_;
    my @domains = $self->search_related("biodata_relationship_objects",
					{ is_location => 1 },
					{ order_by => [qw/rank/] }
	);
    return join("+",map { $_->subject->term->name } @domains);
}

=head2 is_equal_to

 Title   : is_equal_to
 Usage   : $arch->is_equal_to()
 Function: true if input represents the same protein domain
           architecture as the object itself
 Returns : boolean
 Args    : biodata object

           or

           reference to array of hashes. Keys should be named
           the same as the names of biodata_relationship and
           biodata columns

=cut

sub is_equal_to {
    my ($self, $other) = @_;
    if (blessed $other && $other->isa(__PACKAGE__)) {
	return $self->_is_equal_to_object($other);
    } elsif (ref $other eq 'ARRAY') {
	return $self->_is_equal_to_arrayref($other);
    }
    return 0;
}

sub _is_equal_to_object {
    my ($self, $other) = @_;

    return 0 if ($self->biosequence->id ne $other->biosequence->id);

    # Load biodata_relationships
    my @br1 = (
	$self->search_related("biodata_relationship_objects",  undef, { order_by => [qw/rank/] })->all,
	$self->search_related("biodata_relationship_subjects", undef, { order_by => [qw/rank/] })->all
	);
    my @br2 = (
	$other->search_related("biodata_relationship_objects",  undef, { order_by => [qw/rank/] })->all,
	$other->search_related("biodata_relationship_subjects", undef, { order_by => [qw/rank/] })->all
	);
    return 0 if (scalar @br1 != scalar @br2);

    # Compare domains
    for (my $i=0; $i<=$#br1; $i++) {
	my ($subj1, $subj2) = ($br1[$i]->subject,$br2[$i]->subject);
	return 0 if ($subj1->term->name ne $subj2->term->name);
	foreach my $prop (qw(start_coord end_coord strand score)) {
	    my $p1 = $br1[$i]->get_column($prop);
	    my $p2 = $br2[$i]->get_column($prop);
	    return 0 if (!defined $p1 &&  defined $p2); 
	    return 0 if (defined  $p1 && !defined $p2);
	    return 0 if (defined $p1 && defined $p2 && $p1 ne $p2);
	}
    }

    return 1;
}

sub _is_equal_to_arrayref {
    my ($self, $other) = @_;

    # Sequence should be the same
    return 0 if ($self->biosequence->id ne $other->[0]{object}{biosequence}->id);

    # Load biodata_relationships
    my @br1 = (
	$self->search_related("biodata_relationship_objects",  undef, { order_by => [qw/rank/] })->all,
	$self->search_related("biodata_relationship_subjects", undef, { order_by => [qw/rank/] })->all
	);
    return 0 if (scalar @br1 != scalar @$other);

    # Domain data that matters: type (term name), coordinates, strand and score
    for (my $i=0; $i<=$#br1; $i++) {
	my ($subj1,$subj2) = ($br1[$i]->subject,$other->[$i]{subject});
	return 0 if ($subj1->term->name ne $subj2->{term}->name);
	foreach my $prop (qw(start_coord end_coord strand score)) {
	    my $p1 = $br1[$i]->get_column($prop);
	    my $p2 = $other->[$i]{$prop};
	    return 0 if (!defined $p1 &&  defined $p2); 
	    return 0 if ( defined $p1 && !defined $p2);
	    return 0 if ( defined $p1 &&  defined $p2 && $p1 ne $p2);
	}
    }

    return 1;
}

# Finalize Moose
1;
