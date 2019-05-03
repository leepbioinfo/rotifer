# POD documentation - main docs before the code

=head1 NAME

Rotifer::DBIC::AnnotationDB::Role::BioseqPropsRole - 
  accessors to biosequence properties 

=head1 SYNOPSIS

  # Implementing a new parser

  with 'Rotifer::DBIC::AnnotationDB::Role::BioseqPropsRole';

=head1 DESCRIPTION

This role defines the properties (default values, triggers, etc)
for modules that will annotate or access the annotation on raw
sequences (Biosequence table).

It includes a Term object for the sequence type (dna,
protein or rna) and an accessor for the whole sequence
ontology.

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::DBIC::AnnotationDB::Role::BioseqPropsRole;

use strict;
use warnings;
#use autodie qw(:all);
use Carp::Clan qr/^Rotifer::DBIC::AnnotationDB/;
use Data::UUID;
use Moose::Role;

=head2 ATTRIBUTES / ACCESSORS

=head2 sequence_type

 Usage   : Rotifer::DBIC::AnnotationDB->new(sequence_type => "dna")
 Function: set sequence type (dna, rna or polypeptide)
 Returns : Rotifer::DBIC::AnnotationDB::Result::Term
 Value   : string or Term object
 Builder : _default_sequence_type
 Trigger : _sequence_type_trigger
 Default : protein

=cut

has 'sequence_type' => (
    is      => "rw", 
    isa     => "Str|Object",
    lazy    => 1,
    builder => '_default_sequence_type',
    );

sub _default_sequence_type {
    my $self = shift;
    my $type = $self->schema->resultset('Term')->find_or_create
	({ 
	    name        => "polypeptide",
	    definition  => "A polypeptide chain of any length",
	    ontology    => $self->sequence_ontology,
	    is_obsolete => 0,
	 });
    $type->find_or_create_related("term_synonyms",{ synonym => "protein" });
    return $type;
}

# Turn strings into Terms (sequence_types)
sub _sequence_type_trigger {
    my ($self, $new, $old) = @_;
    if (!blessed $new) {
	my $trs = $self->schema->resultset('Term');
	my $term = $trs->search_rs
	    ({ 
		ontology_id => $self->sequence_ontology->id,
		'-or' => {
		    name    => $new,
		    synonym => $new,
		}
	     },
	     { join => "term_synonyms" }
	    );
	if (my $count = $term->count) {
	    if ($count == 1) {
		$new = $term->single;
	    } else {
		croak "Ontology ",$self->sequence_ontology->name," is inconsistent: term $new is ambiguous";
	    }
	} else {
	    $new = $trs->create({ name => $new, ontology => $self->sequence_ontology, is_obsolete => 0 });
	}
	$self->sequence_type($new);
    } else {
	croak "Expected Rotifer::AnnotationDB::Result::Term, got a ",ref($new)
	    unless $new->isa("Rotifer::DBIC::AnnotationDB::Result::Term");
    }
}

=head2 sequence_ontology

 Usage   : $parser->sequence_ontology
 Function: default controlled vocabulary for annotations
 Returns : array reference
 Builder : _default_sequence_ontology
 Trigger : _sequence_ontology_trigger
 Default : sequence

=cut

has 'sequence_ontology' => (
    is       => "rw", 
    isa      => "Str|Object",
    lazy     => 1,
    builder  => '_default_sequence_ontology',
    trigger  => \&_sequence_ontology_trigger,
    );

sub _default_sequence_ontology {
    my $self = shift;
    return $self->schema->resultset('Ontology')->
	find_or_create({ name => "sequence" });
}

sub _sequence_ontology_trigger {
    my ($self, $new, $old) = @_;
    my $default = $self->_default_parser_ontology;
    unless (grep { $_ eq $default } @$new) {
        push(@$new, $default);
        $self->sequence_ontology($new);
    }
}

=head2 METHODS

This role requires/provides no methods.

=cut

1;
