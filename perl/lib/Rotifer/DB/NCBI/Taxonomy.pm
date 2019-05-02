# $Id: Taxonomy.pm,v 1.11 2010/07/14 12:46:52 rfsouza Exp $
#
# BioPerl module for Rotifer::DB::NCBI::Taxonomy
#
# Cared for by Robson Francisco de Souza <rfsouza-at-gmail-dot-com>
#
# Copyright Robson Francisco de Souza
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Rotifer::DB::NCBI::Taxonomy - a tool to query the NCBI taxonomy database

=head1 SYNOPSIS

  # Simplest use case: import gi2taxonomy subroutine

  use Rotifer::DB::NCBI::Taxonomy qw(gi2taxonomy);
  foreach my $tax (gi2taxonomy(@ARGV)) {
      print join("\n",$tax->{gi},$tax->{name}),"\n";
  }

  #....

  # More flexible/quicker option: as an object

  use Rotifer::DB::NCBI::Taxonomy;

  my $util = Rotifer::DB::NCBI::Taxonomy->new(-batch_size => 4000,
                                           -retry => 5);
  $taxutil->submit(@list);
  while (my $data = $taxutil->next_result) {
   # $data is-a reference to a hash (try: keys %$data)
   ...do something...
  }

=head1 DESCRIPTION

This module implements some subroutines to retrieve and process
taxonomy information from Entrez Taxonomy Database.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  http://bugzilla.bioperl.org/

=head1 AUTHOR - Robson Francisco de Souza

Email rfsouza-at-gmail-dot-com

=head1 CONTRIBUTORS

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Rotifer::DB::NCBI::Taxonomy;

use Rotifer;
use Rotifer::DB::Base;
use Rotifer::DB::NCBI qw(id2update);
use XML::Simple;
use Bio::DB::Taxonomy;
use Bio::DB::EUtilities;
use Bio::Tools::EUtilities;
use Bio::Tree::Tree;
use Bio::Root::Root;
use DBI;
use Exporter;
use strict;
use warnings;
use base qw(Exporter Bio::Root::Root);

# Exported methods/tags
our @EXPORT = ();
our @EXPORT_OK = qw(gi2taxonomy);
our %EXPORT_TAGS = ();

# Default values
our $BATCH_SIZE = 7000;
our $RETRY      = 5;

=head2 new

 Title   : new
 Usage   : $tu = Rotifer::DB::NCBI::Taxonomy->new()
 Function: create a Rotifer::DB::NCBI::Taxonomy instance
 Returns : Rotifer::DB::NCBI::Taxonomy object
 Args    : optional arguments are

 -batch_size => # of sequences to retrieve per attempt
                (default 7000)

 -database   => list of NCBI database names to search
                default: protein, nucleotide

 -retry      => # of attempts to reconnect and retrieve
                (default 5)

 -sql        => SQL Database containing gi2taxon table

 -taxdump    => path to directory containing NCBI's 
                Taxonomy database dump files (names.dmp,
                nodes.dmp)

 -fakegis    => file name or reference to an array of
                file names of tables with taxonomy data
                and abbreviations for fakeGI organisms 
                (i.e. organisms whose data does not came
                 from the NCBI)

 -preferred  => reference to an array of lower case taxon
                names to select for short lineage 
                descriptions. This data is used to define
                data that is stored under the "preferred"
                key of the hashes returned by gi2taxonomy
                and next_result 

 -update     => use Rotifer::DB::NCBI::id2update to find the
                newest equivalent GIs to each input ID.
                Enabled by default, set this flag to false
                (0) to disable this feature.

=cut

sub new { 
    my ($class, @args) = shift;
    my $self = $class->SUPER::new(@args);

    # Set defaults and parse options
    my ($bsize, $database, $retry, $sql, $taxdump, $fakeGIs, $update, $preferred) =
	$self->_rearrange([ qw(BATCH_SIZE DATABASE RETRY SQL TAXDUMP FAKEGIS UPDATE PREFERRED) ], @_);
    $self->batch_size(defined $bsize ? $bsize : $BATCH_SIZE );
    $self->retry(defined $retry ? $retry : $RETRY );
    $self->update(defined $update ? $update : 1);
    $self->database(defined $database ? $database : [ "protein", "nucleotide" ]);

    if (defined $sql && ref($sql) eq 'HASH') {
	my ($dsn,$user,$pass) = map { exists $sql->{$_} ? $sql->{$_} : undef } qw(dsn user pass);
	$dsn = "DBI:mysql:database=".$sql->{'database'} if (!defined $dsn && exists $sql->{'database'});
	my $dbh  = DBI->connect($dsn, $user, $pass, { AutoCommit => 0, PrintError => 1, RaiseError => 1}) or die $DBI::errstr;
	$self->{'_dbh'} = $dbh;
    }

    if (defined $taxdump) {
	my $db = Bio::DB::Taxonomy->new(-source    => 'flatfile',
					-directory => $taxdump,
					-nodesfile => "$taxdump/nodes.dmp",
					-namesfile => "$taxdump/names.dmp");
	$self->{'_taxdump'} = $db;
    }

    # Load or set preferred taxons and fakegi organisms
#    $fakeGIs   = "$ENV{ROTIFER_DATA}/taxonomy/fakegi_names.txt" unless (defined $fakeGIs);
    $fakeGIs   = Rotifer::DB::Base->data_path("taxonomy","fakegi_names.txt") unless (defined $fakeGIs);
    unless (defined $preferred) {
	my $path = Rotifer::DB::Base->data_path("taxonomy","taxonomy.txt");
	open(PREF,"<$path") || die "Could not load preferred taxons from $path";
	$preferred = [ map { s/^\s+//; s/\s+$//; $_ } <PREF> ];
	close(PREF);
    }
    $self->_load_fakeGI_organisms(ref $fakeGIs ? @$fakeGIs : $fakeGIs);
    $self->_set_preferred(ref $preferred ? @$preferred : $preferred);

    $self->{'_ids'}     = [];
    $self->{'_results'} = [];
    $self->{'_start'}   = 0;
    return $self;
}

=head2 PREPARE QUERY

=head2 batch_size

 Title   : batch_size
 Usage   : $taxutil->batch_size(1000)
 Function: Get/Set number of IDs to retrieve per connection
 Returns : integer
 Args    : integer
           
=cut

sub batch_size {
    my ($self, $bsize) = @_;
    $self->{'_batch_size'} = $bsize if (defined $bsize);
    $self->{'_batch_size'} = $BATCH_SIZE if (!defined $self->{'_batch_size'});
    return $self->{'_batch_size'};
}

=head2 database

 Title   : database
 Usage   : $taxutil->database([ "protein" ])
 Function: Get/Set number of IDs to retrieve per connection
 Returns : array reference
 Args    : array reference
           
=cut

sub database {
    my ($self, $ref) = @_;
    $self->{'_database'} = $ref if (defined $ref);
    $self->{'_database'} = [ "protein", "nucleotide" ] if (!defined $self->{'_database'});
    return $self->{'_database'};
}

=head2 retry

 Title   : retry
 Usage   : $taxutil->retry(5)
 Function: Get/Set number of attempts to try to retrieve any
           missing data
 Returns : integer
 Args    : integer
           
=cut

sub retry {
    my ($self, $retry) = @_;
    $self->{'_retry'} = $retry if (defined $retry);
    $self->{'_retry'} = $RETRY if (!defined $self->{'_retry'});
    return $self->{'_retry'};
}

=head2 update

 Title   : update
 Usage   : $taxutil->update(1)
 Function: activate/deactivate updating GIs/accession before
           retrieving any taxonomic information
 Returns : boolean
 Args    : boolean
           
=cut

sub update {
    my ($self, $update) = @_;
    $self->{'_update'} = $update if (defined $update);
    return $self->{'_update'};
}

=head2 submit

 Title   : submit
 Usage   : $taxutil->submit(@gis)
 Function: Get/Set the list of identifiers that will be
           used in the next query
 Returns : (integer) number of proteins IDs added
 Args    : (integer) GI numbers

=cut

sub submit {
    my ($self,@ids) = @_;
    if (scalar @ids) {
	warn "New list of IDs submitted while downloading: previous unprocessed queries will be erased!" if (scalar @{ $self->{'_results'} });
	@{ $self->{'_ids'} } = @ids;
	$self->{'_start'} = 0;
    }
    return @{ $self->{'_ids'} };
}

=head2 USER PARAMETERS

=head2 preferred

 Title   : preferred
 Usage   : $taxutil->preferred($lineage)
 Function: Retrieve the preferred name for a taxon
 Returns : string
 Args    : string representing a lineage

 Note    : returns a zero length string if lineage
           does not contain any taxa in the 
           -preferred list given to TaxUtils->new()

=cut

sub preferred {
    my ($self, $lineage) = @_;
#    print join("\n",$self->{'_preferred_taxons'},@{$preferred}),"\n";
    my $preferred = join(">",grep { $_ = lc($_); exists $self->{'_preferred_taxons'}{$_} } grep { !/^celullar/ } split("; ",$lineage));
    return defined $preferred && length $preferred ? $preferred : "";
}

=head2 RETRIEVE DATA

=head2 next_result

 Title   : next_result
 Usage   : $taxutil->next_result()
 Function: Retrieve taxonomy data for the next query
 Returns : hash reference or undef
 Args    : none

=cut

sub next_result {
    my $self = shift;

    my $res = $self->{'_results'};
    if (!scalar @$res) { # Nothing in cache...
	my $start = $self->{'_start'};
	my $ids   = $self->{'_ids'};
	if ($start <= $#{$ids}) { # Something to search for?
	    my $end = $start + $self->batch_size - 1;
	    $end = $end > $#{$ids} ? $#{$ids} : $end;
	    my @input = @{$ids}[$start..$end];
	    @{$res} = $self->_retrieve_taxonomy(@input);
	    $self->{'_start'} += $self->batch_size;
	} else {
	    return ();
	}
    }

    return shift(@{$res});
}

=head2 ID handlers

=head2 gi2taxonomy

 Title   : gi2taxonomy
 Usage   : $taxutil->gi2taxonomy(@gis)
 Function: Get taxonomy for GIs
 Returns : array of hash references
 Args    : (integer) GI numbers

=cut

sub gi2taxonomy {
    my $self = shift if (ref($_[0]) eq __PACKAGE__ || $_[0] eq __PACKAGE__);
    $self = __PACKAGE__->new if (!defined $self || $self eq __PACKAGE__);
    $self->submit(@_);
    my @res = ();
    while (my $hash = $self->next_result) {
	push(@res,$hash);
    }
    return @res;
}

=head2 gi2taxid

 Title   : gi2taxid
 Usage   : $taxutil->gi2taxid(@gis)
 Function: Retrieve taxon IDs for proteins
 Returns : hash of hash references
 Args    : (integer) GI numbers

=cut

sub gi2taxid {
    my ($self, @ids) = @_;

    return () unless (scalar @ids);

    my %data  = (); # Return this hash
    if (defined $self->{'_dbh'}) {
	%data = $self->_gi_from_sql(@ids);
	my @missing = grep { !exists $data{$_} } @ids;
	my %missing = ();
	eval { %missing = $self->_gi_from_eutils(@missing) };
	print STDERR "Problem while trying to determine the taxon ID of each GI (Input: $#ids):\n$@" if ($@);
	%data = (%data, %missing) if (scalar @missing);
    } else {
	eval { %data = $self->_gi_from_eutils(@ids) };
	print STDERR "Problem while trying to determine the taxon ID of each GI (Input: $#ids):\n$@" if ($@);
    }

    return values %data;
}

=head2 INTERNAL METHODS

=head2 _retrieve_taxonomy

 Title   : _retrieve_taxonomy
 Usage   : $taxutil->_retrieve_taxonomy(@gis)
 Function: Retrieve taxon for NCBI's GIs or accessions
 Returns : array of hash references
 Args    : number of times to attempt retrieval (integer)
           GI numbers (array of integers)

=cut

sub _retrieve_taxonomy {
    my ($self, @ids) = @_;

    my @fake = ();
    if (exists $self->{'_fakeGIs_orgs'}) {
	my $fakeGIsre = $self->{'_fakeGIs_orgs_re'};
	my %fake = map { ($_,1) } grep { /^($fakeGIsre)\d+$/o } @ids;
	@ids = grep { !exists $fake{$_} } @ids;
	foreach my $fake ($self->_taxon_for_fakegis(keys %fake)) {
	    if ($fake->{'complete'}) {
		push(@fake, $fake);
		delete $fake{$fake->{'gi'}};
	    } else { # Not found: maybe accession?
		push(@ids, $fake->{'gi'});
	    }
	}
    }

    @ids = $self->_taxon_for_ncbi_ids(@ids);
    return (@ids, @fake);
}

=head2 _taxon_for_ncbi_ids

 Title   : _taxon_for_ncbi_ids
 Usage   : $taxutil->_taxon_for_ncbi_ids(@gis)
 Function: Retrieve taxon for NCBI's GIs or accessions
 Returns : array of hash references
 Args    : GI numbers (array of integers)

=cut

sub _taxon_for_ncbi_ids {
    my ($self, @ids) = @_;

    return () unless (scalar @ids);

    # Attempt to convert input to updated GIs 
    my %updated = ();
    if ($self->update) {
	eval { %updated = id2update(@ids); };
	print STDERR "ERROR: could not apply id2update to input:\n$@" if ($@);
    } else {
	%updated = map { ($_,$_) } @ids;
    }

    my %cache = ();
    my @updated = values %updated;
    my @start = @updated;
    for (my $i=0; $i<=$self->retry; $i++) {
	my @output = $self->gi2taxid(@updated);
	if (exists $self->{'_taxdump'} && $i == 0) {
	    @output = $self->_taxon_from_dump(@output);
	} else {
#	    if (!scalar @output && $#updated < 1000) { # Some fucker crashed EUtils
#		eval { @output = $self->_taxon_from_genbank(@updated) };
#		print STDERR "Attempt: $i IDs: $#start Processing: $#updated Send: $#output Failed to download sequences from GenBank:\n$@" if ($@);
#	    }

	    if (scalar @output) {
		my @missing = grep { !$_->{'complete'} } @output;
		if (scalar @missing) {
		    @output = grep {  $_->{'complete'} } @output;
		    my (@dead, @eutils);
		    map { $_->{'dead'} ? push(@dead, $_) : push(@eutils, $_) } @missing;
		    print STDERR "Attempt: $i IDs: $#start Processing: $#updated Complete: $#output Incomplete: $#missing EUtils: $#eutils GenBank: $#dead\n" if ($self->verbose);
		    eval { @dead = $self->_taxon_from_genbank(map { $_->{'gi'} } @dead) };
		    print STDERR "Attempt: $i IDs: $#start Processing: $#updated Complete: $#output Incomplete: $#missing EUtils: $#eutils GenBank: $#dead Error from Genbank:\n$@" if ($@);

		    if (scalar @eutils) {
			eval { $self->_taxon_from_eutils(@eutils) };
			print STDERR "Attempt: $i IDs: $#start Processing: $#updated Complete: $#output Incomplete: $#missing EUtils: $#eutils GenBank: $#dead Error from eutils:\n$@" if ($@);
			my @missing = grep { !$_->{'complete'} } @eutils;
			if (scalar @missing && $#missing < 1000) {
			    @eutils = grep {  $_->{'complete'} } @eutils;
			    eval { @missing = $self->_taxon_from_genbank(map { $_->{'gi'} } @missing) };
			    print STDERR "Attempt: $i IDs: $#start Processing: $#updated Complete: $#output Incomplete: $#missing EUtils: $#eutils GenBank: $#dead Genbank error :\n$@" if ($@);
			    push(@eutils, @missing);
			}
		    }

		    push(@output, @dead, @eutils);
		    @missing = grep { !$_->{'complete'} } @output;
		    print STDERR "Attempt: $i IDs: $#start Processing: $#updated Complete: $#output Incomplete: $#missing EUtils: $#eutils GenBank: $#dead\n" if ($self->verbose);
		}
	    }
	}

	@output = grep { $_->{'complete'} } @output;
	%cache = (%cache, map { ($_->{'gi'}, $_) } @output);
	@updated = grep { !exists $cache{$_} } @start;
	last if (!scalar @updated);
	print STDERR join(" ","Attempt: $i (failed) IDs: $#start Complete: $#output Missing:",scalar(@updated),"Missing:",@updated),"\n";
    }

    if (scalar @updated) {
	for (my $i=0; $i<=$self->retry; $i++) {
	    my @output = ();
	    eval { @output = $self->_taxon_from_genbank(@updated) };
	    print STDERR "Error while trying to download GenBank file for missing entries:\n$@" if ($@);
	    %cache = (%cache, map { ($_->{'gi'},$_) } @output);
	    @updated = grep { !exists $cache{$_} } @start;
	    last if (!scalar @updated);
	}
    }

    # Restore input order and ids
    my @res = (); my @nogi = (); my @failed = ();
    for (my $i=0; $i<=$#ids; $i++) {
	if (!exists $updated{$ids[$i]}) {
	    push(@nogi, $ids[$i]);
	    next;
	}
	my $id = $updated{$ids[$i]};
	if (!exists $cache{$id}) {
	    push(@failed, $ids[$i]);
	    next;
	}
	my %hash = %{$cache{$id}};
	if ($ids[$i] =~ /^\d+$/) {
	    $hash{'gi'} = $ids[$i];
	} else {
	    $hash{'accession'} = $ids[$i];
	}
	push(@res, { %hash });
    }

    # Last warnings
    if (scalar @failed) {
	print STDERR "!!WARNING!! Taxonomic information could not be found for the following ".scalar(@failed)." queries: ";
	print STDERR join(" ",@failed)."\n";
    }
    if (scalar @nogi) {
	print STDERR "!!WARNING!! The following ".scalar(@nogi)." queries could not be converted to GIs and may not be from the NCBI: ";
	print STDERR join(" ",@nogi)."\n";
    }

    return @res;
}

=head2 _taxon_for_fakegis

 Title   : _taxon_for_fakegis
 Usage   : $taxutil->_taxon_for_fakegis(@gis)
 Function: Retrieve taxons for fakeGIs
 Returns : hash of hash references
 Args    : (integer) GI numbers

=cut

sub _taxon_for_fakegis {
    my ($self, @ids) = @_;

    my $orgs = $self->{'_fakeGIs_orgs'};
    my @taxonomy = ();
    foreach my $key (@ids) {
	(my $organism = $key) =~ s/\d+$//; # here is were this subroutine parses fakeGIs
	if (exists $orgs->{$organism}) {
	    push(@taxonomy,
		 {
		     gi        => $key,
		     accession => $key,
		     abbreviation => $organism,
		     dead      => 0,
		     taxid     => undef,
		     lineage   => $orgs->{$organism}{'lineage'},
		     preferred => $orgs->{$organism}{'lineage'},
		     name      => $orgs->{$organism}{'name'},
		     complete  => 1,
		 });
	} else {
	    push(@taxonomy,
		 {
		     gi        => $key,
		     accession => $key,
		     abbreviation => undef,
		     dead      => 0,
		     taxid     => 0,
		     lineage   => 'NO_LINEAGE',
		     preferred => 'NO_LINEAGE',
		     name      => 'NO_NAME',
		     complete  => 0,
		 });
	}
    }

    return @taxonomy;
}

=head2 _taxon_from_eutils

 Title   : _taxon_from_eutils
 Usage   : $taxutil->_taxon_from_eutils(@gis)
 Function: Actually retrieve data from TaxonomyDB
 Returns : hash of hash references
 Args    : (integer) GI numbers

=cut

sub _taxon_from_eutils {
    my ($self, @data) = @_;

    return () unless (scalar @data);

    my %taxon = (0 => { 'lineage' => 'NO_LINEAGE', 'name' => 'NO_NAME' });
    map { $taxon{$_->{'taxid'}} = { %{$taxon{"0"}} } } @data;

    # Epost TaxIds: send list
    my $factory = Bio::DB::EUtilities->new(-eutil      => 'epost',
					   -email      => 'nobody@nowhere.com',
					   -db         => 'taxonomy',
					   -id         =>  [ grep { $_ } keys %taxon ],
#					   -verbose    => 1,
					   -usehistory => 'y');
    my $parser = Bio::Tools::EUtilities->new(-eutil => 'epost', -response => $factory->get_Response);
    my $history = $parser->next_History;
    return () if (!defined $history);

    # Efetch for TaxIds: retrieve XML descriptions
    $factory = Bio::DB::EUtilities->new(-eutil   => 'efetch',
					-email   => 'nobody@nowhere.com',
					-db      => 'taxonomy',
					-report  => 'xml',
#					-verbose => 1,
					-history => $history);
    my $xml = XMLin($factory->get_Response->content);
    my @taxa = ref($xml->{'Taxon'}) eq 'ARRAY' ? @{$xml->{'Taxon'}} : ref($xml) eq 'HASH' ? $xml->{'Taxon'} : die "Unable to parse taxonomy entry: unknown response";
    foreach my $taxon (@taxa) {
	my $taxid = $taxon->{'TaxId'} || 0;
	$taxon{$taxid}->{'lineage'} = $taxon->{'Lineage'}        || 'NO_LINEAGE';
	$taxon{$taxid}->{'name'}    = $taxon->{'ScientificName'} || 'NO_NAME';
    }

    # Collect final results
    foreach my $data (@data) {
	my $taxid = $data->{'taxid'};
	if (exists $taxon{$taxid}) {
	    $data->{'name'}      = $taxon{$taxid}->{'name'};
	    $data->{'lineage'}   = $taxon{$taxid}->{'lineage'};
	    $data->{'preferred'} = $self->preferred($data->{'lineage'});
	    $data->{'complete'}  = 1;
	}
    }

    return @data;
}

=head2 _taxon_from_genbank

 Title   : _taxon_from_genbank
 Usage   : $taxutil->_taxon_from_genbank(@gis)
 Function: Retrieve taxon by downloading and
           parsing Genbank files
 Returns : hash of hash references
 Args    : (integer) GI numbers

=cut

sub _taxon_from_genbank {
    my ($self, @ids) = @_;

    return () unless (scalar @ids);

    use File::Which qw(which);
    use File::Temp qw(tempfile);
    my ($fh, $filename) = tempfile(UNLINK => 1);
    print $fh join("\n",@ids),"\n";
    close($fh);
    my $edirect = which("edirect.pl") || "edirect.pl";
    open(FETCH,"cat $filename | xargs -n1 $edirect -fetch -format gp -db protein -id |") || do {
	warn "Failed to run edirect.pl";
	return ();
    };

    use Bio::SeqIO;
    my $io = Bio::SeqIO->new(-fh     => \*FETCH,
			     -format => "genbank");
    my %found = ();
    my @data = ();
    while (my $seq = $io->next_seq) {
	my $spc = $seq->species;
	if (defined $spc) {
	    push(@data,{ gi        => $seq->primary_id,
			 accession => $seq->accession_number || "",
			 abbreviation => undef,
			 dead      => 0,
			 taxid     => $spc->ncbi_taxid || 0, 
			 lineage   => join("; ",reverse $spc->classification) || "NO_LINEAGE",
			 name      => $spc->binomial || "NO_NAME",
			 complete  => 1,
		 });
	    $data[$#data]->{'lineage'}   = "cellular organisms; " . $data[$#data]->{'lineage'} if
		($data[$#data]->{lineage} !~ /^cellular/ && $data[$#data]->{lineage} =~ /^(bacteria|archaea|eukaryote)/i);
	    $data[$#data]->{'preferred'} = $self->preferred($data[$#data]->{'lineage'});
	} else {
	    print STDERR "Could not determine species data for GI ".$seq->primary_id;
	    push(@data,{ gi        => $seq->primary_id,
			 accession => $seq->accession_number, 
			 abbreviation => undef,
			 dead      => 0,
			 taxid     => 0,
			 lineage   => "NO_LINEAGE",
			 preferred => "NO_LINEAGE",
			 name      => "NO_NAME",
			 complete  => 1,
		 });
	}
	$found{$data[$#data]->{'gi'}} = 1;
    }
    $io->close;
    close(FETCH);

    if (my @missing = grep { !exists $found{$_} } @ids) {
	print STDERR join(" ","The following IDs were not found in GenBank:",@missing),"\n" if ($self->verbose);
    }

    return @data;
}

=head2 _taxon_from_dump

 Title   : _taxon_from_dump
 Usage   : $taxutil->_taxon_from_dump(@gis)
 Function: Retrieve taxon data from NCBI's
           TaxonomyDB dumps
 Returns : hash of hash references
 Args    : (integer) GI numbers

=cut

sub _taxon_from_dump {
    my ($self, @data) = @_;
    my $db = $self->{'_taxdump'};

    return () unless (scalar @data);

    my %taxon = map { ($_->{'taxid'},0) } @data;
    my $tree_functions = Bio::Tree::Tree->new();
    foreach my $taxid (keys %taxon) {
	my $tax = $db->get_taxon($taxid);
	my @lineage = $tree_functions->get_lineage_nodes($tax);
	$taxon{$taxid} = join("; ",map { $_->node_name } grep { $_->node_name !~ /^cellular/ } @lineage);
    }
    foreach my $data (@data) {
	$data->{'lineage'}   = $taxon{$data->{'taxid'}} || 0;
	$data->{'preferred'} = $self->preferred($data->{'lineage'});
    }

    return @data;
}

=head2 _gi_from_eutils

 Title   : _gi_from_eutils
 Usage   : $taxutil->_gi_from_eutils(@gis)
 Function: Find the GI of each accession using EUtils
 Returns : hash of hash references
 Args    : (integer) GI numbers

=cut

sub _gi_from_eutils {
    my ($self, @ids) = @_;
    my %data  = (); # Return this hash

    return () unless (scalar @ids);

    foreach my $db (@{$self->database}) {
	# EPost GIs: send list
	my $factory = Bio::DB::EUtilities->new(-eutil      => 'epost',
					       -db         => $db,
					       -id         =>  \@ids,
					       -email      => 'nobody@nowhere.com',
#					       -verbose    => 1,
					       -usehistory => 'y');
	my $parser = Bio::Tools::EUtilities->new(-eutil => 'epost', -response => $factory->get_Response);
	my $history = $parser->next_History;

	# ESummary for GIs: retrieve docsums
	$factory = Bio::DB::EUtilities->new(-eutil   => 'esummary',
					    -email   => 'nobody@nowhere.com',
					    -db      => $db,
#					    -verbose => 1,
					    -history => $history);
	$parser = Bio::Tools::EUtilities->new(-eutil => 'esummary',
#					      -verbose => 1,
					      -response => $factory->get_Response);

	# Parse docsums
	while (my $docsum = $parser->next_DocSum) {
	    my ($gi) = $docsum->get_contents_by_name('Gi');
	    my ($accession) = $docsum->get_contents_by_name('Caption');
	    if (($docsum->get_contents_by_name('Status'))[0] eq 'dead') {
		#warn "Sequence entry ".$docsum->get_id." is dead!";
		$data{$gi} = {
		    gi        => $gi,
		    accession => $accession,
		    abbreviation => undef,
		    taxid     => 0,
		    dead      => 1,
		    lineage   => 'NO_LINEAGE',
		    preferred => 'NO_LINEAGE',
		    name      => 'NO_NAME',
		    complete  => 0,
		};
	    } else {
		my ($taxid) = $docsum->get_contents_by_name('TaxId');
		$data{$gi} = { 
		    gi        => $gi,
		    accession => $accession,
		    abbreviation => undef,
		    taxid     => $taxid,
		    dead      => 0,
		    lineage   => 'NO_LINEAGE',
		    preferred => 'NO_LINEAGE',
		    name      => 'NO_NAME',
		    complete  => 0,
		};
	    }
	}

	# Deleting incomplete entries
	my @missing = ();
	foreach my $id (@ids) {
	    if (exists $data{$id}) {
		delete $data{$id} if (!defined $data{$id}->{'taxid'});
	    } else {
		push(@missing, $id);
	    }
	}

	# Try another database if something is missing
	if (scalar @missing) {
	    @ids = @missing;
	} else {
	    last;
	}
    } # foreach my $db (@{$self->database})

    return %data;
}

=head2 _gi_from_sql

 Title   : _gi_from_sql
 Usage   : $taxutil->_gi_from_sql(@gis)
 Function: Map accessions to GIs using a gi2taxon
           table in a SQL server
 Returns : hash of hash references
 Args    : (integer) GI numbers

=cut

sub _gi_from_sql {
    my ($self, @ids) = @_;
    return () unless (scalar @ids);

    my $where = join(" OR ", map { "gi = '$_'" } @ids);
    my $sth   = $self->{'_dbh'}->prepare_cached("SELECT gi, ncbi_taxon_id as taxid FROM gi2taxon WHERE $where");
    $sth->execute();
    my %data = ();
    while (my $row = $sth->fetchrow_hashref) {
	$data{$row->{'gi'}} = {
	    gi        => $row->{'gi'},
	    accession => 'unknown',
	    abbreviation => undef,
	    taxid     => $row->{'taxid'},
	    dead      => 0,
	    lineage   => 'NO_LINEAGE',
	    preferred => 'NO_LINEAGE',
	    name      => 'NO_NAME',
	    complete  => 0,
	};
    }

    # Retrieve taxon data from NCBI using e-utils
    eval { $self->_taxon_from_eutils(values %data) };

    return %data;
}

=head2 _load_fakeGI_organisms

 Title   : _load_fakeGI_organisms
 Usage   : $taxutil->_load_fakeGI_organisms(@gis)
 Function: Load table of organisms whose sequences
           are identified by fake GIs
 Returns : array of hash references
 Args    : one or more table file names

=cut

sub _load_fakeGI_organisms {
    my ($self, @files) = @_;

    # Parse fakegi organisms data
    my %orgs = ();
    my @columns = qw(organism file name preferred_taxons);
    foreach my $file (@files) {
	open(my $fh,"<$file");
	while (<$fh>) {
	    chomp;
	    my @data = map { s/^\s+//; s/\s+$//; $_ } split(/\s+:\s+/);
	    $orgs{$data[0]} = {}; # Create empty anonymous hash
	    map { $orgs{$data[0]}->{$columns[$_]} = $data[$_] } (1..3);
	    $orgs{$data[0]}->{'lineage'} = $data[3];
	}
	close($fh);
    }

    # Store regular expression and table
    $self->{'_fakeGIs_orgs'} = { %orgs };
    $self->{'_fakeGIs_orgs_re'} = join("|",keys %{$self->{'_fakeGIs_orgs'}});

    return %orgs;
}

=head2 _fastacmd

 Title   : _fastacmd
 Usage   : $taxutil->_fastacmd(@gis)
 Function: 
 Returns : 
 Args    : 

=cut

sub _fastacmd {
    my $self = shift;
    my %ids  = ();
    my @ids  = map { $ids{$_} = 1; $_ } @_;

    return () unless (scalar @ids);

    use File::Temp qw(tempfile);
    my ($fh, $file) = tempfile();
    print $fh join("\n",@ids),"\n";
    close($fh);

    open(FASTACMD,"fastacmd -d nr -T T -i $file |");
    my $gi = undef;
    my %data = ();
    while (<FASTACMD>) {
	chomp;
	/^NCBI sequence id: (\S+)/ && do {
	    my @a = split(/\|/, $1);
	    $gi = $a[1];
	    unless (exists $ids{$gi}) {
		$gi = undef;
		next;
	    }
	    $a[3] =~ s/\.\d+$//;
	    $data{$gi} = { gi => $gi, accession => $a[3] };
	    next;
	};
	next unless (defined $gi);

	/^NCBI taxonomy id: (\S+)/ && do {
	    $data{$gi}->{'taxid'} = $1;
	    next;
	};

	/Scientific name: (.+)/ && do {
	    $data{$gi}->{'name'} = $1;
	    next;
	};
    }
    close(FASTACMD);
    unlink($file) if ( -f $file );

    return %data;
}

=head2 _set_preferred

 Title   : _set_preferred
 Usage   : $taxutil->_set_preferred()
 Function: Store list of preferred taxon names
 Returns : 
 Args    : ARRAY reference

=cut

sub _set_preferred {
    my ($self, @preferred) = @_;
    $self->{'_preferred_taxons'} = { map { (lc($_),1) } @preferred };
}

sub DESTROY {
    my ($self) = shift;
    $self->{'_dbh'}->disconnect if (exists $self->{'_dbh'});
}

1;
