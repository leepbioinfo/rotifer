# POD documentation - main docs before the code

=head1 NAME

Rotifer::DB - high level Rotifer methods to access databases

=head1 SYNOPSIS

  # Retrieving FASTA sequences
  use Rotifer::DB qw(id2fasta);
  my $fasta = id2fasta($gi);

=head1 DESCRIPTION

Rotifer::DB aggregates methods to retrieve
data from both remote and local databases. 

=head1 CONFIGURATION

Rotifer::DB depends on configuration parameters that are set
in a file called "db.yml" that is searched for the user's
~/.rotifer/etc directory or, if not found, in the etc/rotifer
directory under your Rotifer's installation path.

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::DB;

use strict;
use warnings;
use autodie qw(:all);
use Carp::Clan qr/^Rotifer::DB/;
use File::Monitored qw(:temp);
use List::MoreUtils qw(uniq);
use Rotifer::Utils qw(nr2ids);
use Rotifer::DB::NCBI;
use Rotifer::DB::FADB;
use Scalar::Util qw(reftype);
use Sub::Exporter -setup => {
    exports => [qw/id2fasta/],
};

=head2 id2fasta

 Title   : id2fasta
 Usage   : Rotifer::DB->id2fasta(\%opts, @gis)
 Function: get sequences or sequence data 
 Returns : filename or file handle
 Args    : array of sequence identifiers (strings)

=cut

sub id2fasta {
    my $self = shift if (defined $_[0] && ($_[0] eq __PACKAGE__ || UNIVERSAL::isa($_[0],__PACKAGE__)));
    my $opts = { %{shift(@_)} } if (defined $_[0] && UNIVERSAL::isa($_[0],"HASH"));
    return undef unless (scalar @_);
    my %ids = map { ($_,0) } @_; # Register queries

    # Output
    my ($outfh,$output);
    my $output_isa_file = 1;
    if (exists $opts->{output} && defined $opts->{output}) {
	$output = delete $opts->{output};
	my $reftype = reftype($output);
	if (defined $reftype && $reftype eq "GLOB") {
	    $outfh = $output;
	    $output = undef;
	    $output_isa_file = 0;
	} else { # A filename was given
	    open($outfh,">$output");
	}
    } else { # Create temporary file
	($outfh, $output) = &ttempfile();
    }

    # Choose sequences source
    if (defined $opts->{"db"}) {
	my $flag = 0;
	my $indexfile = $opts->{"db"} . ".idx";
	my @blastdbs = glob($opts->{"db"}."\.{p,n}??");
	if ( -f $opts->{"db"} && !scalar(@blastdbs)) {
	    # Take care of the index file
	    use Bio::Index::Fasta;
	    use Bio::SeqIO;
	    use File::stat;
	    carp join(" ","Retrieving",scalar(@_),"sequences from local (indexed) FASTA file",$opts->{"db"}) if ($opts->{debug});

	    # Check pre-existing index
	    if ( -f $indexfile ) {
		my $indexstat = stat($indexfile);
		my $filestat  = stat($opts->{"db"});
		if ($indexstat->mtime < $filestat->mtime) {
		    unlink($indexfile);
		    $flag = 1;
		}
	    } else {
		$flag = 1;
	    }

	    # Build/load index
	    my $index = eval { Bio::Index::Fasta->new(-filename=>$indexfile, -write_flag=>$flag) };
	    if ($@) {
		warn $@ if ($opts->{debug} > 3);
		$flag = 1;
		unlink($indexfile);
		$index = Bio::Index::Fasta->new(-filename=>$indexfile, -write_flag=>$flag);
	    }
	    if ($flag) {
		carp join(" ","Indexing local FASTA file",$opts->{"db"}) if ($opts->{debug});
		$index->id_parser(\&_get_id);
		$index->make_index($opts->{"db"});
	    }

	    # Retrieve and print
	    my $out = Bio::SeqIO->new(-fh=>$outfh, -format=>"fasta");
	    foreach my $id (@_) {
		my $seqobj = $index->fetch($id);
		next unless (defined $seqobj);
		$out->width($seqobj->length) if ($opts->{clean});
		$out->write_seq($seqobj);
		$ids{$id} = 1;
	    }
	} # if (defined $opts->{"db"}  && !scalar(@blastdbs))
    } # if (defined $opts->{"db"})

    # Isolate fake GIs and retrieve from other/external databases
    my @fakeGI = ();
    my @other = ();
    foreach my $id (grep { !$ids{$_} } @_) {
	if (Rotifer::DB::FADB->is_fakeGI($id)) {
	    push(@fakeGI,$id);
	} else {
	    push(@other,$id);
	}
    }
    Rotifer::DB::NCBI->id2fasta({ %$opts, output => $outfh }, @other);
    Rotifer::DB::FADB->fakegi2fasta({ %$opts, output => $outfh }, @fakeGI);

    # Merge results and return
    close($outfh) if ($output_isa_file);
    return $output;
}

sub _get_id {
    my $header = shift;
    use Rotifer::Utils qw(nr2ids);
    my (@ids) = map { $_->{accession} } nr2ids($header);
    return (@ids);
}

# Make perl compiler happy
1;
