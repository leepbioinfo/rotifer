=head1 NAME

Rotifer::DBIC::AnnotationDB::ResultSet::Biodata - biodata resultset

=head1 DESCRIPTION
 
This modules extends DBIx::Class::ResultSet for the table biodata.

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::DBIC::AnnotationDB::ResultSet::Biodata;

use strict;
use warnings;
use Carp::Clan qr/^Rotifer::DBIC/;
use Moose;
use MooseX::NonMoose;
use MooseX::MarkAsMethods autoclean => 1;
use Rotifer::DB::NCBI::Taxonomy;
use Rotifer::Utils;
use Scalar::Util qw(blessed reftype);
extends 'DBIx::Class::ResultSet';

=head2 ATTRIBUTES / ACCESSORS

=head1 Utilities

=head2 retrieve_all_organisms

 Title   : retrieve_all_organisms
 Usage   : $count = $self->retrieve_all_organisms(\%opts,@biodata)
 Function: downlaod NCBI taxonomy information for biodata
 Returns : array of Rotifer::DBIC::AnnotationDB::ResultSet::Organism
 Args    : options for Rotifer::DB::NCBI::Taxonomy (hash reference)
           and

           array references of arguments to be passed to
           DBIx::Class::ResultSet::search(). These arguments will be
           used to query for biodata that needs updating.

           and/or

           array of Rotifer::DBIC::AnnotationDB::Result::Biodata

=cut

sub retrieve_all_organisms {
    my ($self, $opts, @data) = @_;

    # Some pre-requisites (prefer newest GIs)
    my $ors = $self->result_source->schema->resultset("Organism");

    # Query biodata and get best ids
    my $count = 0;
    my %id2biodata = ();
    foreach my $data (@data) {
	next unless (defined $data);
	my @bd = ($data);
	@bd = $self->search(@$data)->all if (!blessed $data && reftype($data) eq "ARRAY");
	foreach my $biodata (@bd) {
	    next unless $biodata->isa("Rotifer::DBIC::AnnotationDB::Result::Biodata");
	    next if (!$opts->{"-update"} && $biodata->has_column_loaded("organism_id") && defined $biodata->get_column("organism_id"));
	    my @xrefs = sort {
		my $aacc = $a->accession;
		my $bacc = $b->accession;
		$Rotifer::Utils::SeqIdTypeRank{$a->acctype} <=> $Rotifer::Utils::SeqIdTypeRank{$b->acctype} ||
		    (looks_like_number($aacc) && looks_like_number($bacc) ? $aacc <=> $bacc : $aacc cmp $bacc) ||
		    $b->version <=> $a->version
	    } $biodata->dbxrefs->search({ acctype => { '-in' => [ keys %Rotifer::Utils::SeqIdTypeRank ] } });
	    next unless (scalar @xrefs);
	    next if (exists $id2biodata{ $xrefs[0]->accession });
	    $id2biodata{ $xrefs[0]->accession } = $biodata;
	    $count++;
	}
    }

    # Retrieve taxonomy
    carp "Update organisms for $count biodata: retrieving NCBI's taxonomic data...\n" if ($self->result_source->schema->debug);
    my $tutil = Rotifer::DB::NCBI::Taxonomy->new(%{$opts});
    $tutil->submit(keys %id2biodata);
    my %orgs = ();
    while (my $hashref = $tutil->next_result) {
	my $biodata = 
	    exists $id2biodata{$hashref->{gi}}        ? $id2biodata{$hashref->{gi}}        :
	    exists $id2biodata{$hashref->{accession}} ? $id2biodata{$hashref->{accession}} :
	    undef;
	next unless (defined $biodata);
	$hashref->{description} = "missing taxon data" if ($hashref->{name} eq 'NO_NAME');

	# Organism
	my $orgID = 
	    exists $hashref->{taxid}        && defined $hashref->{taxid}        ? $hashref->{taxid}        :
	    exists $hashref->{abbreviation} && defined $hashref->{abbreviation} ? $hashref->{abbreviation} :
	    do { carp "Can't insert organism without NCBI Taxon ID and/or abbreviation"; next };

	my $organism = $orgs{$orgID};
	if (!defined $organism) {
	    $organism = $ors->update_or_create
		({
		    ncbi_taxon_id  => $hashref->{taxid},
		    abbreviation   => $hashref->{abbreviation},
		    name           => $hashref->{name},
		    lineage        => $hashref->{preferred},
		    classification => $hashref->{lineage},
		    comment        => $hashref->{comment} || undef,
		 },
		 {
		     key => exists $hashref->{taxid} && defined $hashref->{taxid} ?
			 'organism_ncbitaxonid_uniq' :
			 'organism_abbreviation_uniq'
		 });
	    $orgs{$orgID} = $organism;
	}

	$biodata->organism($organism);
	$biodata->update;
    } # while (my $hashref = $tutil->next_result)

    my @orgs = sort { $a->id <=> $b->id } values %orgs;
    return @orgs;
}

__PACKAGE__->meta->make_immutable();
1;

