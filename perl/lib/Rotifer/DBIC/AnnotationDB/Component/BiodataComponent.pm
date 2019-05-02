=head1 NAME

Rotifer::DBIC::AnnotationDB::Component::BiodataComponent - methods for biosequences

=head1 DESCRIPTION
 
Extensions to Rotifer::DBIC::Result::Biodata.

=head1 Class inheritance

This module inherits or consumes the following classes/interfaces/roles:

  Rotifer::DBIC::AnnotationDB::Component::WithDataset
  DBIx::Class
  Moose::Obejct
  Bio::SeqI
  Bio::Root::Root

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::DBIC::AnnotationDB::Component::BiodataComponent;

use strict;
use warnings;
use Carp::Clan qr/^Rotifer::DBIC::AnnotationDB/;
use Data::UUID;
use DBIx::Class::ResultClass::HashRefInflator;
use Moose;
use Memoize;
use Rotifer::Utils qw(ids2nr);
use Rotifer::DB::NCBI::Taxonomy;
use Scalar::Util qw(blessed);
extends qw(
 Rotifer::DBIC::AnnotationDB::Component::WithDataset
 DBIx::Class
 Moose::Object
 Bio::SeqI
 Bio::Root::Root
);
#memoize('display_id');

=head1 Utility methods

=head2 retrieve_organism

 Title   : retrieve_organism
 Usage   : $organism = $self->retrieve_organism
 Function: retrieve and reset taxonomy information
 Returns : 
 Args    : Rotifer::DB::NCBI::Taxonomy options
 Note    : this method downloads data from NCBI

=cut

sub retrieve_organism {
    my ($self, %taxopts) = @_;
    carp "Update organism information for biodata: downloading taxonomic data from NCBI...\n";

    # Some pre-requisites
    my $ors = $self->result_source->schema->resultset('Organism');

    # Use newest GIs or accessions for each biodata in the database
    my @acc = $self->search_related("Dbxref",
				    { acctype => "GI" },
				    { columns => [ "accession" ] }
	)->all;
    return unless (scalar @acc);

    # Retrieve taxonomy: use only the newest GI's organism data
    my $tutil = Rotifer::DB::NCBI::Taxonomy->new(%taxopts);
    my ($hashref) = sort { $b->{gi} <=> $a->{gi} } $tutil->gi2taxonomy(@acc);
    return unless (defined $hashref);
    $hashref->{comment} = "missing taxon data" if ($hashref->{name} eq 'NO_NAME');
    my $organism = $ors->find_or_new({ name           => $hashref->{name},
				       defined $hashref->{abbreviation} ? (abbreviation   => $hashref->{abbreviation}) : (),
				       defined $hashref->{taxid}        ? (ncbi_taxon_id  => $hashref->{taxid})        : (),
				       lineage        => $hashref->{preferred},
				       classification => $hashref->{lineage},
				       comment        => $hashref->{comment} || undef,
				     }, { key => 'organism_uniq' });
    if ($organism->in_storage) { # Update name
	$organism->name($hashref->{name});
	$organism->update;
    } else {
	$organism->insert;
    }

    # Add taxonomy link
    $self->organism($organism);
    $self->update;

    return 1;
}

# =head1 Bio::SeqFeatureI-compatible methods

=head2 location

 Title   : location
 Usage   : $rs = $biodata->location()
 Function: get all locations for this biodata
 Returns : DBIx::Class::ResultSet
 Args    : none

=cut

sub location {
    my $self = shift;

    # Check if this biodata was registered as the base of some coordinate system
    #
    # Coordinate systems are biodata entries that either represent a sequence (is_rawsequence == 1)
    # or have a biodata_relationship to themselves (subject_id = object_id) and is_location = 1
    my $location = $self->search_related_rs("biodata_relationship_subjects",
					    { subject_id       => { "=" => \'object_id' },
					      'predicate.name' => 'identity',
					      is_location      => 1,
	                                    },
			      	            { join     => 'predicate',
                                              order_by => [qw/start_coord end_coord rank/],
                                            }
       ); #'}});
    return $location if ($location->count);

    # Retrieve a simple location
    $location = $self->search_related_rs("biodata_relationship_subjects",
					 { subject_id  => { "!=" => \'object_id' },
				           is_location => 1,
	                                 },
       		  	                 {
				           order_by => [ qw/start_coord end_coord rank/ ],
                                         },
	); #'}});

    return $location;
}

=head1 Bio::FeatureHolderI methods

=head2 add_SeqFeature

 Title   : add_SeqFeature
 Usage   : $rs->add_SeqFeature($biodata, 'EXPAND', $relationship);
 Function: add sub-feature (related biodata)
 Returns : feature added
 Args    : up to three arguments
           1 - (obligatory) Biodata row object
           2 - 'EXPAND' or undef: adjust the parent coordinates 
           3 - (optional) hash reference with values to other
                          BiodataRelationship columns

=cut

sub add_SeqFeature {
    my ($self, $biodata, $expand, $brel) = @_;
    my $brs  = $self->schema->resultset("Biodata");
    my $brrs = $self->schema->resultset("BiodataRelationship");
    my $ors  = $self->schema->resultset("Ontology");
    my $trs  = $self->schema->resultset("Term");

    # Set predicate for relationship
    my $predicate = $brel->{predicate} || "part_of";
    if (!ref $predicate) {
	$predicate = {
	    name        => $predicate, 
	    ontology    => $ors->find_or_create({ name => "relationships" }),
	    is_obsolete => 0,
	};
    }
    $brel->{predicate} = $trs->find_or_create($predicate) if (!blessed $predicate);

    # Autocreate child and set is_location property
    $brel->{subject} =  blessed $biodata ? $biodata : $brs->find_or_create($biodata);
    $brel->{is_location} = 0 unless (exists $brel->{is_location});

    # Set rank
    if (!exists $brel->{rank}) {
	my $brel->{rank} = $self->search_related_rs("biodata_relationship_objects", $brel)->get_column("rank")->max;
    }

    # Add relationship
    $brel = $self->find_or_create_related("biodata_relationship_objects", $brel);

    # Update limits
    carp "Support EXPAND not implemented" if (defined $expand);

    return $brel;
}

=head2 get_SeqFeatures

 Title   : get_SeqFeatures
 Usage   : my @feats = $seq->get_SeqFeatures();
 Function: retrieve just the toplevel sequence features attached to this seq
 Returns : array of biodata
 Args    : none

=cut

sub get_SeqFeatures {
    my ($self) = shift;
    my $rs = $self
	->search_related_rs("biodata_relationship_objects")
	->search_related_rs("subject");
    if (defined $self->get_column("biosequence_id")) {
	my $rs2 = $self->result_source->schema->resultset("Biodata")
	    ->search_rs({ is_obsolete    => 0,
			  is_rawsequence => 1,
			  biosequence_id => $self->get_column("biosequence_id"),
			})
	    ->search_related_rs("biodata_relationship_objects")
	    ->search_related_rs("subject");
	$rs->union($rs2);
    }
    return $rs->all;
}

=head2 get_all_SeqFeatures

 Title   : get_all_SeqFeatures
 Usage   : @features = $annseq->get_all_SeqFeatures()
 Function: returns all SeqFeatures, included sub SeqFeatures
 Returns : an array of Bio::SeqFeatureI objects
 Args    : none

=cut

sub get_all_SeqFeatures {
    my $self = shift;
    my @all  = ();
    my @feat = $self->get_SeqFeatures;
    my %done = map { ($_->id,$_) } @feat;
    while (scalar @feat) {
	push(@all, @feat);
	@feat = grep { !$done{$_->id} } map { $_->get_SeqFeatures } @feat;
	map { $done{$_->id} = $_ } @feat;
    }
    return @all;
}

=head2 feature_count

 Title   : feature_count
 Usage   : $seq->feature_count()
 Function: Return the number of SeqFeatures attached to a sequence
 Returns : integer representing the number of SeqFeatures
 Args    : none

=cut

sub feature_count {
    return scalar(shift->get_all_SeqFeatures);
    return shift->search_related_rs("biodata_relationship_objects", {}, { columns => [ qw/object_id/ ], distinct => 1 })->count;
}

=head1 Bio::SeqI-compatible methods

=head2 primary_seq

 Title   : primary_seq
 Usage   : $obj->primary_seq($newval)
 Function: Retrieve the underlying Bio::PrimarySeqI object if available.
           This is in the event one has a sequence with lots of features
           but want to be able to narrow the object to just one with
           the basics of a sequence (no features or annotations).
 Returns : Bio::PrimarySeqI
 Args    : Bio::PrimarySeqI or none;

       See Bio::PrimarySeqI for more information

=cut

{ no warnings 'once';
  *primary_seq = \&Rotifer::DBIC::AnnotationDB::Result::Biodata::biosequence;
}

=head1 Bio::PrimarySeqI-compatible methods

=cut

=head2 display_id, primary_id, accession_number

 Mostly useless but implemented!

=cut

sub display_id {
    my $self = shift;
    my $drs = $self->dbxrefs({ acctype => { "-in" => \@Rotifer::Utils::SeqIdTypes } }, { cache => 1 });
    return $self->get_column("identifier") unless ($drs->count);
    my $group_rs = $drs->search({},{ columns => [qw/accgroup/], distinct => 1, cache => 1 });
    $group_rs = $group_rs->get_column("accgroup")->max_rs->as_query;
    my $d = $self->dbxrefs->search(undef, { accgroup => $group_rs, cache => 1 });
    $d->result_class("DBIx::Class::ResultClass::HashRefInflator");
    my @ids = $d->all;
    map { delete $_->{dbxref_id}; $_->{description} = "" } @ids;
    return $self->get_column("identifier") unless (scalar @ids);
    my $display_id = ids2nr({ description => 0 }, @ids);
    return $display_id;
}

{ no warnings 'once';
#  *display_id       = \&Rotifer::DBIC::AnnotationDB::Result::Biodata::identifier;
  *display_name     = *display_id;
  *accession        = \&Rotifer::DBIC::AnnotationDB::Result::Biodata::identifier;
  *accession_number = \&Rotifer::DBIC::AnnotationDB::Result::Biodata::identifier;
}

=head2 seq

 Title   : seq
 Usage   : $seq = $rs->seq("acgt")
 Function: get residues for the related biosequence
 Returns : string (scalar)
 Args    : none

=cut

sub seq {
    my ($self) = @_;
    my $bioseq = $self->biosequence;
    return defined $self->biosequence ? $self->biosequence->residues : undef;
}

=head2 subseq

 Title   : subseq
 Usage   : $text = $rs->subseq(10,20)
 Function: Returns the subseq from start to end, where the first base
           is 1 and the number is inclusive, i.e. 1-2 are the first two
           bases of the sequence.

           Start cannot be larger than end but can be equal.

 Returns : string or undef
 Args    : pair of integers

=cut

sub subseq {
    my $seq = shift->biosequence;
    return defined $seq ? $seq->subseq(@_) : undef;
}

=head2 trunc

 Title   : trunc
 Usage   : $text = $rs->trunc(10,20)
 Function: same as subseq() but return a new object
 Returns : Rotifer::DBIC::AnnotationDB::Result::Biodata
 Args    : pair of integers

 Note    : only DBIx::Class::Row->new() is executed
           you must run $obj->insert() to store it

=cut

sub trunc {
    my ($self) = shift;
    return undef if ($self->length == 0 || !defined $self->biosequence);
    my %data = $self->get_inflated_columns;
    return $self->result_source->schema->resultset("Biodata")->find_or_new({ %data,
									     biosequence => $data{biosequence}->subseq(@_),
									     identifer   => "$data{identifier}:" . join("..",@_),
									   });
}

=head2 length

 Title   : length
 Usage   : $text = $rs->length()
 Function: get the sequence length
 Returns : integer or undef
 Args    : none

=cut

sub length {
    my $self = shift;
    my $length = 0;
    my $loc = $self->location;
    my $count = $loc->count;
    if ($count) {
	my @obj = $loc->search({},{ columns => [qw/ object_id /], distinct => 1 })->all;
	if (scalar(@obj) > 1) {
	    carp "Length is not defined for biodata projected against multiple coordinate systems";
	    return 0;
	} else {
	    my $lastend = 0;
	    foreach my $l ($loc->all) {
		my $start = $l->start_coord;
		$start = $lastend + 1 if ($start < $lastend + 1); # Correct for overlaps
		next if ($l->end_coord < $start);
		$length += $l->end_coord - $start + 1;
		$lastend = $l->end_coord;
	    }
	}
    }
    if (defined $self->get_column("biosequence_id")) {
	my $seqlen = $self->biosequence->length;
	if ($length > $seqlen) {
	    carp "Data inconsistency: biodata length > sequence length for ",$self->identifier;
	    return undef;
	}
	$length = $self->biosequence->length;
    }
    return $length;
}

=head2 description

 Title   : description
 Usage   : $text = $rs->description()
 Function: does nothing
 Returns : undef
 Args    : 

=cut

{ no warnings 'once';
  *desc = \&Rotifer::DBIC::AnnotationDB::Result::Biodata::description;
};

=head2 alphabet

 Title   : alphabet
 Usage   : $text = $rs->alphabet()
 Function: does nothing
 Returns : undef
 Args    : 

=cut

sub alphabet {
    my $seq = shift->biosequence;
    return defined $seq ? $seq->alphabet : undef;
}

=head2 is_circular

 Title   : is_circular
 Usage   : $text = $rs->is_circular()
 Function: find whether this sequence is a circular replicon
 Returns : boolean
 Args    : none

=cut

sub is_circular {
    my $seq = shift->biosequence;
    return defined $seq ? $seq->is_circular : undef;
}

__PACKAGE__->meta->make_immutable(inline_constructor => 0);
1;
