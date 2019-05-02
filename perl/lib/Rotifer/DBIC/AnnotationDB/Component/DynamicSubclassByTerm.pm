=head1 NAME

Rotifer::DBIC::AnnotationDB::Component::DynamicSubclassByTerm - 
 auto-detect subclasses for tables typed by term

=head1 DESCRIPTION

This DBIx:component implements dynamic subclassing of row objects for tables whose
types are defined by a term_id column that refers to entries in the term table.

Based on DBIx::Class::DynamicSubclass and the dynamic subclassing examples in
DBIx::Class::Cookbook.

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::DBIC::AnnotationDB::Component::DynamicSubclassByTerm;

use Carp::Clan qr/^Rotifer::DBIC::AnnotationDB/;
use List::MoreUtils qw(uniq);
use base qw(DBIx::Class);
use strict;
use warnings;

our $VERSION = '0.01000';

=head2 classify

 Title   : classify
 Usage   : $rs->classify
 Function: redefine object's class based on its type (term)
 Returns : reblessed object (also in-place)
 Args    : none

 Note    : this method will do nothing unless there is a subclass of
           Biodata whose name matches the name or a synonym of the
           referred Term table entry. One example:

           Rotifer::DBIC::AnnotationDB::Result::Biodata::reference

           is implemented and Biodata entries of type "reference" are
           properly reblessed.

=cut

sub classify {
    my $self = shift;

    # In case there is some subclass
    if ($self->has_column_loaded("term_id") && defined $self->term) {
	my $class_name = ref $self ? ref $self : $self;
	my $term_name  = $self->term->name;
	if (my @invalid = ($term_name =~ /([^\w\_\-\s])/g)) {
	    @invalid = uniq @invalid;
	    if ($self->result_source->schema->can("debug") && $self->result_source->schema->debug > 1) {
		carp "Term name $term_name contains the invalid characters: @invalid.\n";
		carp "Please use only alphanumeric characters, underscores, dashes ('-') and/or spaces for names.";
	    }
	    return $self;
	}
	$term_name =~ s/^\s+//;
	$term_name =~ s/\s+$//;
	$term_name =~ s/\s/_/g;
	$term_name =~ s/\-+/_/g;

	my $start = 1;
	my @syn = ("${class_name}::${term_name}");
	while (my $class = shift @syn) {
	    if ($self->ensure_class_found($class)) {
		$self->ensure_class_loaded($class);
		bless $self, $class;
		return $self;
	    }
	    if ($start) { # Look at synonyms too...
		@syn = sort grep { !/[^\w:]/ } map {
		    my $synonym = $_->synonym;
		    $synonym =~ s/^\s+//;
		    $synonym =~ s/\s+$//;
		    $synonym =~ s/\-/_/g;
		    $synonym =~ s/\s/_/g;
		    "${class_name}::${synonym}";
		} $self->term->term_synonyms;
		$start = 0;
	    }
	}
    }

    return $self;
}

=head1 OVERLOADED METHODS

These section documents DBIx::Class::Row methods that are overloaded to
implement automatic reblessing of Biodata objects based on the term_id
column.

=head2 new

 Title   : new
 Usage   : $rs->new
 Function: new row objects (see DBIx::Class::Row)
 Returns : row object
 Args    : same as DBIx::Class::Row::new

=cut

sub new {
    my $this = shift;
    my $data = shift;

    my $deferred;
    if ($this->can('add_frozen_columns')) {
        my $real_columns = $this->result_source_instance->_columns;
        map {
            $deferred->{$_} = delete $data->{$_}
                unless index($_, '-') == 0 or exists $real_columns->{$_};
        } keys %$data;
    }

    my $ret = $this->next::method($data, @_);
    $ret->classify;

    if ($deferred) {
        $ret->set_columns($deferred);
    }

    return $ret;
}

=head2 inflate_result

 Title   : inflate_result
 Usage   : $rs->inflate_result
 Function: convert data to row objects (see DBIx::Class::Row)
 Returns : row object
 Args    : same as DBIx::Class::Row::inflate_result

=cut

sub inflate_result {
    my $self = shift;
    my $ret = $self->next::method(@_);
    $ret->classify;
    return $ret;
}

sub store_column {
    my ($self, $column, $value) = @_;

    if ($column eq "term_id") {
        my $ret = $self->next::method($column, $value);
        $self->classify;
        return $ret;
    }

    $self->next::method($column, $value);
}

1;
