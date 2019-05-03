# POD documentation - main docs before the code

=head1 NAME

Rotifer::DB::NCBI::EUtils - simplyfied interface for Bio::DB::EUtilities
                         and Bio::Tools::EUtilities

=head1 SYNOPSIS

  # Convert protein identifiers to DNA GI's
  use Rotifer::DB::NCBI::EUtils qw(elink);
  my @table = elink({ from => "protein", to => "nucleotide" }, $gi);

=head1 DESCRIPTION

Rotifer::DB::NCBI::EUtils provides a simpler interfaces to retrieve
data from NCBI E-Utils service on top of the interface implemented
in Chris J. Fields's Bio::DB/Tools::EUtilities.

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 DEPENDENCIES

=over

=item Bio::DB::EUtilities

=item Bio::Tools::EUtilities

=back

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::DB::NCBI::EUtils;

use strict;
use warnings;
use autodie qw(:all);
use Bio::DB::EUtilities;
use Bio::Tools::EUtilities;
use Carp::Clan qr/^Rotifer::DB::NCBI::EUtils/;
use File::Monitored qw(:temp);
use IO::String;
#use IO::File;
use List::MoreUtils qw(minmax);
use LWP::UserAgent;
use Scalar::Util qw(blessed reftype);
use Rotifer::Utils qw(nr2ids);
use XML::Simple;
use Sub::Exporter -setup => {
    exports => [qw/id2genbank id2fasta
                   get_seqids accession2gi gi2accession
                   elink efetch epost esummary/],
};

# Global variables
our $MAX_DOWNLOAD_ATTEMPTS = 10;
our @NON_BIOPERL_OPTIONS   = qw/count debug epost ouput retry sleep/;
our $DEFAULT_SERVER        = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils';
our %DEFAULT_SERVICES      = ( 'elink' => 'elink.fcgi' );

=head2 id2fasta

 Title   : id2fasta
 Usage   : Rotifer::DB::NCBI::EUtils->id2fasta(\%opts, @gis)
 Function: retrieve sequences in FASTA format
 Returns : file name or file handle (same as key "output")
 Args    : array of sequence identifiers (strings)

 Users may precede the list of identifiers with an
 optional anonymous hash reference that may include
 any key supported by Bio::DB::EUtilities::new and

   debug => print progress messages for debugging
  output => file name or file handle
   retry => numbers of attempts to download the data

=cut

sub id2fasta {
    my $self = shift if (defined $_[0] && ($_[0] eq __PACKAGE__ || ref($_[0]) eq __PACKAGE__));
    my $opts = { %{shift(@_)} } if (defined $_[0] && ref($_[0]) eq "HASH");
    my @ids  = @_;

    # Process options
    $opts->{rettype} = 'fasta';
    $opts->{retmode} = 'text';
    my $attempts = $opts->{retry} || $MAX_DOWNLOAD_ATTEMPTS;
    my $debug    = exists $opts->{debug} ? $opts->{debug} : 0;

    # Ouput
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

    # Prepare ids
    my %ids = ();
    my %gi = map { @{$_} } get_seqids(grep { /^\d+$/ } @ids);
    for (my $i=0; $i<=$#ids; $i++) {
	$ids[$i] = $gi{$ids[$i]} if (exists $gi{$ids[$i]});
	$ids{$ids[$i]} = 1;
    }

    # Query, retrieve, process...
    my $round  = 0;
    my $last = "";
    while (scalar(@ids) && $round < $attempts) {
	carp join(" ","id2fasta -> Server error, missing",scalar(@ids),". Redo #$round:",@ids) if ($debug > 1 && $round > 0);
	carp join(" ","id2fasta -> Retrieving",scalar(@ids),"of",scalar(@_),"in round #",$round) if ($debug);
	my $filename = efetch($opts, @ids);
#	++$round && next unless (defined $filename);
	last unless (defined $filename);

	open(FETCHED, "<$filename");
	my $load = 0;
	while (<FETCHED>) {
	    next if /^\s*$/;
	    if (/^>/) {
		$load = 0;
		my @hash = grep { exists $ids{$_->{accession}} } nr2ids($_);
		if (scalar @hash) {
		    map {
			delete $ids{$_->{accession}};
			delete $ids{$_->{accession} .".". $_->{version}};
		    } @hash;
		    $load = 1;
		}
		$_ = "\n$_" if (length $last && $last !~ /\n$/);
	    } else {
		s/\s//g if ($opts->{clean});
	    }
	    print $outfh $_ if ($load);
	    $last = $_;
	}
	print $outfh "\n" if ($load && length $last && $last !~ /\n$/);
	close(FETCHED);

	$round++;
	@ids = grep { exists $ids{$_} } @ids;
    }
    close($outfh) if ($output_isa_file);

    # Complain about things that did not work
    if (scalar @ids) {
	carp join(" ","id2fasta -> Server error, missing",scalar(@ids),":",@ids) if ($debug);
    }

    return $output;
}

=head2 id2genbank

 Title   : id2genbank
 Usage   : Rotifer::DB::NCBI::EUtils->id2genbank(\%opts, @gis)
 Function: download genbank files
 Returns : file name or file handle
 Args    : array of sequence identifiers (strings)

 Users may precede the list of identifiers with an
 optional anonymous hash reference that may include
 any key supported by Bio::DB::EUtilities::new plus

   debug => print progress messages for debugging
  output => file name or file handle to dump output
   retry => numbers of attempts to download the data

=cut

sub id2genbank {
    my $self = shift if (defined $_[0] && ($_[0] eq __PACKAGE__ || ref($_[0]) eq __PACKAGE__));
    my $opts = { %{shift(@_)} } if (defined $_[0] && ref($_[0]) eq "HASH");
    my @ids  = @_;

    # Process options
    $opts->{rettype} = 'gbwithparts';
    $opts->{retmode} = 'text';
    my $retry = $opts->{retry} || $MAX_DOWNLOAD_ATTEMPTS;
    my $debug = exists $opts->{debug} ? $opts->{debug} : 0;

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
    #print STDERR join(" ","I2G",$outfh,$output,@ids),"\n";

    # Prepare ids
    my %ids = ();
    my %gi = map { @{$_} } get_seqids(grep { /^\d+$/ } @ids);
    for (my $i=0; $i<=$#ids; $i++) {
	$ids[$i] = $gi{$ids[$i]} if (exists $gi{$ids[$i]});
	$ids{$ids[$i]} = 1;
    }

    # Query, retrieve, process...
    my $round  = 0;
    while (scalar(@ids) && $round < $retry) {
	carp join(" ","id2genbank -> Server error, missing",scalar(@ids),". Redo #$round:",@ids) if ($debug > 1 && $round > 0);
	carp join(" ","id2genbank -> Retrieving",scalar(@ids),"of",scalar(@_),"in round #",$round) if ($debug);
	my $temp = efetch($opts, @ids);
	++$round && next unless (defined $temp && -s $temp);

	# Fetching
	my @found  = ();
	my $string = "";
	open (FETCHED,"<$temp") || die "Could not open temporary file $temp";
	while (<FETCHED>) {
	    next if /^\s*$/;

	    # End of sequence: store
	    if (/^\/\/\s*$/) {
		if (length $string && $string =~ /^LOCUS/) {
		    $string .= $_;
		    print $outfh $string;
		    map { delete $ids{$_} } @found;
		}
		$string = "";
		@found = ();
		next;
	    }

	    # Accession
	    if (/^ACCESSION\s+(\S+)/) {
		push(@found,$1);
	    }

	    # Accession.Version
	    if (/^VERSION\s+(\S+)/) {
		push(@found,$1);
	    }

	    # GI
	    elsif (/^VERSION\s+\S+\s+GI:(\d+)/) {
		push(@found,$1);
	    }

	    $string .= $_; # Concatenate Genbank file
	} # while (<FETCHED>)
	close(FETCHED);
	unlink($temp);

	@ids = grep { exists $ids{$_} } @ids;
	$round++;
    } # while  (scalar(keys %ids) && $round < $retry)
    close($outfh) if ($output_isa_file);

    # Complain about things that did not work
    if (scalar @ids) {
	carp join(" ","id2genbank -> Server error, missing",scalar(@ids),":",@ids) if ($debug);
    }

    return $output;
}

=head2 get_count

 Title   : get_count
 Usage   : get_count( epost(@ids) )
 Function: get the number of results for epost(s)
 Returns : array of integers
 Args    : Bio::DB::EUtilities::History

=cut

sub get_count {
    my $self = shift if (defined $_[0] && ($_[0] eq __PACKAGE__ || ref($_[0]) eq __PACKAGE__));
    my $opts = shift if (defined $_[0] && ref($_[0]) eq "HASH");

    my @count = ();
    foreach my $history (@_) {
	my %hash = (-eutil  => 'esearch',
		    -email  => 'nobody@nowhere.com',
		    _hashref_to_bioperl($opts),
		    -db      => $opts->{db} || 'protein',
		    -retmax  => 0,
		    -history => $history,
	    );
	my $factory = Bio::DB::EUtilities->new(%hash);
 	push(@count,$factory->get_count);
    }

    return @count;
}

=head2 get_seqids

 Title   : get_seqids
 Usage   : Rotifer::DB::NCBI::EUtils->get_seqids(@ids)
 Function: get accession numbers and GIs for any ids
 Returns : array of arrays (GI and accession pairs)
 Args    : array of sequence identifiers (strings)
 Alias   : accession2gi and gi2accession

 Users may precede the list of identifiers with a
 reference to an anonymous hash with the keys

   retry => numbers of attempts to download the data
   debug => print progress messages for debugging

=cut

sub get_seqids {
    my $self = shift if (defined $_[0] && ($_[0] eq __PACKAGE__ || ref($_[0]) eq __PACKAGE__));
    my $opts = (defined $_[0] && ref($_[0]) eq "HASH") ? shift : {};
    #my @ids = map { s/\.\d+$// if /^[^\.]+\.\d+$/; $_ } @_; 
    my @ids = @_;

    # Parse options: copy $opts to avoid clashes with other subroutines 
    my $retry = $opts->{retry} || $MAX_DOWNLOAD_ATTEMPTS;
    my $debug = exists $opts->{debug} ? $opts->{debug} : 0;
    $opts = {
	%$opts,
	rettype => "seqid",
	db      => 'protein',
	retmode => 'text',
    };

    # Absolutely necessary!
    delete $opts->{output};
    #$opts->{post} = 0;

    # Send/retrieve
    my @array = ();
    my $round = 0;
    my (%gi,%acc) = ();
    map { if (/^\d+$/) { $gi{$_} = 1 } else { $acc{$_} = 1 } } @ids;
    while (scalar(@ids) && $round < $retry) { # Shit happens... try again!
	carp join(" ","get_seqids -> missing",scalar(@ids),". Redo #$round:",@ids) if ($debug > 1 && $round > 0);
	carp join(" ","get_seqids -> retrieving ",scalar(@ids),"of",scalar(@_),"in round #",$round) if ($debug);
	#my $accfile = scalar(keys %acc) ? _efetch_only($opts, keys %acc) : undef;
	my $accfile = efetch($opts, keys %acc);
	my $gifile  = efetch($opts,keys %gi);

	# Fetching
	my $acc = undef;
	my $version = undef;
	foreach my $file ($gifile, $accfile) {
	    next unless (defined $file);
	    open (FETCHED, "<$file" );
	    while (<FETCHED>) { # Parse ASN.1
		next if /^\s*$/;
		if (/accession\s+\"(\S+)\"/) {
		    $acc = $1;
		} elsif (/version\s+(\d+)/) {
		    $version = $1;
		} elsif (/Seq-id\s+::=\s+gi\s+(\d+)/) {
		    next unless (defined $acc);
		    delete $acc{$acc};
		    $acc .= ".$version" if (defined $version);
		    push(@array, [ $1, $acc ]);
		    delete $acc{$acc};
		    delete $gi{$1};
		    undef($version);
		    undef($acc);
		};
	    }
	    close(FETCHED);
	}

	@ids = grep { exists $gi{$_} || exists $acc{$_} } @ids;
	$round++;
    } # while  (scalar(keys %ids) && $round < $retry)

    # Complain about things that did not work
    if (scalar @ids) {
	carp join(" ","get_seqids -> Server error, missing ids:",scalar(@ids),"ids:",@ids,
		  "\nIt seems at least one of the accessions above is not from any NCBI sequence database\n");
    }

    return @array;
}
{no warnings 'once';
 *accession2gi = \&get_seqids;
 *gi2accession = \&get_seqids;
}

=head2 esummary

 Title   : esummary
 Usage   : Rotifer::DB::NCBI::EUtils->esummary(\%opts, @gis)
 Function: use the esummary service
 Returns : array of Bio::Tools::EUtilities::Summary::DocSum
 Args    : array of sequence identifiers (strings)

 Users may precede the list of identifiers with an
 optional anonymous hash reference that may include
 any key supported by Bio::DB::EUtilities::new

=cut

sub esummary {
    my $self = shift if (defined $_[0] && ($_[0] eq __PACKAGE__ || ref($_[0]) eq __PACKAGE__));
    my $opts = defined ($_[0]) && ref($_[0]) eq "HASH" ? { %{ shift(@_) } } : {};

    # Isolate histories, replace accessions with GIs and drop unknown accessions
    my %gi  = map { ($_->[1],$_->[0]) } get_seqids(grep { !/^\d+$/ } @_);
    my @ids = grep { /^\d+$/ } map { exists $gi{$_[$_]} ? $gi{$_[$_]} : $_[$_] } 0..$#_;
    #my @ids = @_;
    return () if (!scalar @ids);

    # Epost: send GI list
    my @history = epost($opts, @ids);

    # Retrieve summaries
    my @docsum = ();
    for (my $i=0; $i<=$#history; $i++) {
	# Build factory
	my %hash = (-eutil   => 'esummary',
		    -db      => $opts->{db} || 'protein',
		    -email   => 'nobody@nowhere.com',
		    -history => $history[$i],
		    _hashref_to_bioperl($opts));
	my $factory = Bio::DB::EUtilities->new(%hash);

	# Iterate over docsums
	while (my $docsum = $factory->next_DocSum) {
	    push(@docsum, $docsum);
	}
    }

    return @docsum;
}

=head2 epost

 Title   : epost
 Usage   : Rotifer::DB::NCBI::EUtils->epost(\%opts, @gis)
 Function: use the epost service
 Returns : Bio::Tools::EUtilities::
 Args    : array of sequence identifiers (strings)

 Users may precede the list of identifiers with an
 optional anonymous hash reference that may include
 any key supported by Bio::DB::EUtilities::new

=cut

sub epost {
    my $self = shift if (defined $_[0] && ($_[0] eq __PACKAGE__ || ref($_[0]) eq __PACKAGE__));
    my $opts = shift if (defined $_[0] && ref($_[0]) eq "HASH");
    # Epost: send list
    my $factory = Bio::DB::EUtilities->new(-eutil   => 'epost',
					   -email   => 'nobody@nowhere.com',
					   -db      => $opts->{db} || 'protein',
					   -id      =>  [ @_ ],
#					   -verbose => 1,
					   -usehistory => 'y',
					   _hashref_to_bioperl($opts));
    my $response = $factory->get_Response;
    my $parser   = Bio::Tools::EUtilities->new(-eutil => 'epost', -response => $response);
    my $history  = $parser->next_History;
    return $history;
}

=head2 elink

 Title   : elink
 Usage   : Rotifer::DB::NCBI::EUtils->elink(\%opts, @gis)
 Function: use the elink service
 Returns : array of arrays
 Args    : array of sequence identifiers (strings)

 This subroutine returns a table with three columns:

 1 - input identifiers
 2 - output identifier
 3 - comma separated list of link names

 Users may precede the list of identifiers with an
 optional anonymous hash reference with keys

    debug => print progress messages for debugging
     dump => do not parse results, return array of XML
   retmax => maximum number of queries per download

 Any other options will be added to the URL sent to
 EUtils. See http://www.ncbi.nlm.nih.gov/books/NBK25499/
 for details.

 Example:

 elink({ linkname => "protein_structure" }, @ARGV)

 Default values for required ELink options are:

       db => nucleotide
   dbfrom => protein

 Note: not based on Bio::DB::EUtilities
       Uses LWP::Simple and XML::Simple

 See

 http://eutils.ncbi.nlm.nih.gov/entrez/query/static/entrezlinks.html

 for a list of link names

=cut

sub elink {
    my $self = shift if (defined $_[0] && ($_[0] eq __PACKAGE__ || ref($_[0]) eq __PACKAGE__));
    my $opts = shift if (defined $_[0] && ref($_[0]) eq "HASH");
    my %fopts = %$opts; # Copy options, just in case...
    my $ua = LWP::UserAgent->new;

    # Accession to GI
    my %gi  = map { ($_->[1],$_->[0]) } get_seqids(grep { !/^\d+$/ } @_) if (!exists $fopts{'dbfrom'} || $fopts{'dbfrom'} =~ /sequence|protein|^nuc(core|leotide)/);
    my @ids = grep { /^\d+$/ } map { exists $gi{$_[$_]} ? $gi{$_[$_]} : $_[$_] } 0..$#_;
    #my @ids = @_;

    # Process other options
    my $debug   = delete $fopts{"debug"}  || 0;
    my $dump    = delete $fopts{"dump"}   || 0;
    my $retmax  = delete $fopts{"retmax"} || 500;
    my $retry   = delete $fopts{"retry"}  || 5;
    my $sleep   = delete $fopts{"sleep"}  || 0;
    $fopts{dbfrom} = 'protein'    unless (exists $fopts{dbfrom});
    $fopts{db}     = 'nucleotide' unless (exists $fopts{db});
    $fopts{idtype} = 'acc'        unless (exists $fopts{idtype}); # support for accession.version

    # Prepare base URL
    my $base   = "$DEFAULT_SERVER/$DEFAULT_SERVICES{elink}?";
    $base .= join("&",map { $_ . "=" . $fopts{$_} } keys %fopts);

    # Process IDs
    my $attempt = 1;
    my @table   = ();
    my %table   = ();
    my $start   = 0;
    my $end     = $#ids > $retmax-1 ? $retmax-1 : $#ids;
    while ($start <= $#ids) {
	sleep $sleep if ($start > 0); # Wait a few seconds before next download...
	carp join(" ","elink -> retrieving links for chunk ${start}..$end of",scalar(@ids)) if ($debug);

	# Submit chunk and parse
	my @chunk = @ids[$start..$end];
	my $url = $base.join("",map { "&id=$_" } @chunk);
	carp "GET $url" if ($debug > 1);
	my $response = $ua->get($url);
	if (!$response->is_success) {
	    carp "Download attempt $attempt failed! Status:\n" . $response->status_line if ($debug);
	    if ($attempt < $retry) { # lets try again!
		$attempt++;
		next;
	    }
	    $response = undef;
	}

	# You want raw XML?
	if ($dump) {
	    push(@table, $response->content);
	}

	# Process links
	elsif (defined $response) {
	    my $xml = XMLin($response->content);
	    foreach my $linkset (ref $xml->{LinkSet} eq "ARRAY" ? @{$xml->{LinkSet}} : $xml->{LinkSet}) {
		my $from = $linkset->{IdList}{Id};

		# If the target database is not in NCBI, $from == undef
		next unless (defined $from);

		# Process each set
		my $set = ref $linkset->{LinkSetDb} eq "ARRAY" ? $linkset->{LinkSetDb} : [ $linkset->{LinkSetDb} ];
		my %worst = (); my %rejected = ();
		foreach my $link (@$set) {
		    my $type = $link->{LinkName};
		    next unless (defined $type);
		    next unless (exists $link->{Link});

		    # Select
		    my $results = ref($link->{Link}) eq "ARRAY" ? $link->{Link} : [ $link->{Link} ];
		    foreach my $id (@$results) {
			my $to = $id->{Id};
			next if (!defined $to);
			$table{$from}{$to}{$type} = exists $id->{Score} ? $id->{Score} : -1;
		    }
		} # foreach my $link (@$set)
	    } # foreach my $linkset (ref $xml->{LinkSet} eq "ARRAY" ? @{$xml->{LinkSet}} : $xml->{LinkSet})
	} # elsif (defined $response)

	# Prepare for next chunk
	$start = $end+1;
	$end   = $end+$retmax;
	$end   = $#ids if ($#ids < $end);
	$attempt = 1;
    } # while ($start <= $#ids)

    # Flatten hash of hashes
    if (!$dump) {
	foreach my $from (sort keys %table) {
	    my @to = keys %{$table{$from}};
	    next if (!scalar @to);
	    foreach my $to (@to) {
		my @keys = sort keys %{$table{$from}{$to}};
		my $sum = undef;
		my @values = map { $sum += $table{$from}{$to}{$_} if ($table{$from}{$to}{$_} >= 0); $table{$from}{$to}{$_} } @keys;
		push(@table,[ $from, $to, join(", ",@keys) ]);
		push(@{ $table[$#table] }, join(", ",@values)) if (defined $sum);
	    }
	}
    }

    return @table;
}

=head2 efetch

 Title   : efetch
 Usage   : Rotifer::DB::NCBI::EUtils->efetch(\%opts, @gis)
 Function: use the efetch service
 Returns : path to file with all data
 Args    : array of sequence identifiers (strings)
           or Bio::Tools::EUtilities::History objects

 Users may precede the list of identifiers with an
 optional anonymous hash reference that may include
 any key supported by Bio::DB::EUtilities::new plus

   count => anonymous array with history size
   debug => print progress messages for debugging
   epost => boolean: whether to use epost if @_ > 200
  output => file name or file handle to dump output
   retry => numbers of attempts to download the data

=cut

sub efetch {
    my $self = shift if (defined $_[0] && ($_[0] eq __PACKAGE__ || ref($_[0]) eq __PACKAGE__));
    my $opts = shift if (defined $_[0] && ref($_[0]) eq "HASH");
    if ((scalar(@_) < 250 && !grep { blessed $_ } @_) || (exists $opts->{epost} && !$opts->{epost})) {
	return _efetch_only($opts, @_);
    } else {
	return _epost_efetch($opts, @_);
    }
}

=head1 Internal methods

=head2 _epost_efetch

 Title   : _epost_efetch
 Usage   : Rotifer::DB::NCBI::EUtils->_epost_efetch(\%opts, @gis)
 Function: epost then download with efetch
 Returns : path to file with all data
 Args    : array of sequence identifiers (strings)
           or Bio::Tools::EUtilities::History objects

 Users may precede the list of identifiers with an
 optional anonymous hash reference that may include
 any key supported by Bio::DB::EUtilities::new plus

   count => anonymous array with history size
   debug => print progress messages for debugging
  output => file name or file handle to dump output
   retry => numbers of attempts to download the data

=cut

sub _epost_efetch {
    my $self = shift if (defined $_[0] && ($_[0] eq __PACKAGE__ || ref($_[0]) eq __PACKAGE__));
    my $opts = { %{ shift(@_) } } if (defined $_[0] && ref($_[0]) eq "HASH");

    # Isolate histories, replace accessions with GIs and drop unknown accessions
    my @history = grep { blessed $_ } @_;
    my @ids = grep { !blessed $_ } @_;
    my @count = get_count(@history) if (scalar @history);

    # Epost: send GI list
    if (scalar @ids) {
	push(@history, epost($opts, @ids));
	push(@count, $#ids);
    }
    return () if (!scalar @history);

    # Parse options
    my $retry    = $opts->{retry} || $MAX_DOWNLOAD_ATTEMPTS;
    my $retmax   = exists $opts->{retmax}   ? $opts->{retmax}   : 5000;
    my $retstart = exists $opts->{retstart} ? $opts->{retstart} : 0;
    my $debug    = exists $opts->{debug}    ? $opts->{debug}    : 0;

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

    # Iterate over histories
    my $i = 0;
    foreach my $history (@history) {
	# Prepare efetch
	my %hash = (-eutil   => 'efetch',
		    -db      => $opts->{db} || 'protein',
		    -rettype => 'gbwithparts',
		    -retmode => 'text',
		    -email   => 'nobody@nowhere.com',
		    -history => $history,
		    _hashref_to_bioperl($opts));
	my $factory = Bio::DB::EUtilities->new(%hash);

	# Retrieve each chunk
	my $attempt = 0;
      RETRIEVE_SEQS:
	while ($retstart <= $count[$i]) {
	    carp join(" ","_epost_efetch -> Chunk start",$retstart,"with",scalar(@ids),"queries.") if ($debug);
	    $factory->set_parameters(-retmax => $retmax, -retstart => $retstart);

	    eval{ # Safely attempt a download
		$factory->get_Response(-cb => 
				       sub { 
					   my ($data) = @_;
					   if ($data =~ /<\?xml.*>.*<ERROR>(.+)<\/ERROR>/) {
					       $@ = $1;
					   } else {
					       $data =~ s/\n\s*\n/\n/;
					       print $outfh $data;
					   }
				       });
	    };

	    if ($@) { # Recover up to $retry times
		confess "_epost_efetch -> Server error: $@.  Try again later" if $attempt == $retry;
		if ($@ =~ /Bad Request/) {
		    carp "_epost_efetch -> query execution failed for chunk $retstart: ".
			$factory->parameter_base->to_request->uri."\n".
			"_epost_efetch -> query execution failed for chunk $retstart:$@"
			if ($debug);
		    return undef;
		} else {
		    carp "_epost_efetch -> Server error, redo #$attempt chunk $retstart:\n" if ($debug);
		}
		$attempt++ && redo RETRIEVE_SEQS;
	    }

	    $retstart += $retmax;
	} # while ($retstart <= $#ids)

	$i++;
    } # foreach my $history (@history)
    close($outfh) if ($output_isa_file);

    return $output;
}

=head2 _efetch_only

 Title   : _efetch_only
 Usage   : Rotifer::DB::NCBI::EUtils->_efetch_only(\%opts, @gis)
 Function: download with efetch
 Returns : path to file with all data or undef (see note below)
 Args    : array of sequence identifiers (strings)

 Users may precede the list of identifiers with an
 optional anonymous hash reference that may include
 any key supported by Bio::DB::EUtilities::new plus

   count => anonymous array with history size
   debug => print progress messages for debugging
  output => file name or file handle to dump output
   retry => numbers of attempts to download the data

 Note: when the optional hash's key 'output' is set to
       to a glob this subroutine will assume the value
       points to an open filehandle and will attempt to
       print all downloaded content to this filehandle.
       It will also return 'undef' instead of a file path.

=cut

sub _efetch_only {
    my $self = shift if (defined $_[0] && ($_[0] eq __PACKAGE__ || ref($_[0]) eq __PACKAGE__));
    my $opts = { %{ shift(@_) } } if (defined $_[0] && ref($_[0]) eq "HASH");

    # Ignore histories, efetch can accept both accessions and GIs
    my @ids = grep { !blessed $_ } @_;

    # Parse options
    my $retry  = $opts->{retry} || $MAX_DOWNLOAD_ATTEMPTS;
    my $retmax = exists $opts->{retmax} ? $opts->{retmax} : 400;
    my $debug  = exists $opts->{debug}  ? $opts->{debug}  : 0;

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

    # Retrieve each chunk
    my $attempt = 0;
    my $start   = 0;
    my $end     = $retmax-1 > $#ids ? $#ids : $retmax-1;
  RETRIEVE_SEQS:
    while ($start <= $#ids) {
	carp join(" ","_efetch_only -> attempt",$attempt,"Chunk",$start,"..",$end,"with",scalar(@ids),"queries.") if ($debug);

	# Prepare efetch
	my %hash = (-eutil   => 'efetch',
		    -db      => $opts->{db} || 'protein',
		    -rettype => 'gbwithparts',
		    -retmode => 'text',
		    -email   => 'nobody@nowhere.com',
		    -id      => [ @ids[$start..$end] ],
		    _hashref_to_bioperl($opts));
	my $factory = Bio::DB::EUtilities->new(%hash);
	carp "_efetch_only -> query for chunk $start .. $end is\n".$factory->parameter_base->to_request->uri."\n" if ($debug > 1);

	eval{ # Safely attempt a download
	    $factory->get_Response(-cb => 
				   sub { 
				       my ($data) = @_;
				       if ($data =~ /<\?xml.*>.*<ERROR>(.+)<\/ERROR>/) {
					   $@ = $1;
				       } else {
					   $data =~ s/\n\s*\n/\n/;
					   print $outfh $data;
				       }
				   });
	};

	if ($@) { # Recover up to $retry times
	    confess "_efetch_only -> Server error:\n$@\nTry again later..." if $attempt == $retry;
	    if ($@ =~ /Bad Request/) {
		carp "_efetch_only -> query execution failed for chunk $start: ".
		    $factory->parameter_base->to_request->uri."\n".
		    "_epost_efetch -> query execution failed for chunk $start:$@"
		    if ($debug);
		carp "_efetch_only -> Query execution failed:$@" if ($debug);
		return undef;
	    } else {
		carp "_efetch_only -> Server error, redo #$attempt chunk $start:\n" if ($debug);
	    }
	    $attempt++ && redo RETRIEVE_SEQS;
	}

	$start += $retmax;
	$end = $start + $retmax - 1;
	$end = $#ids if ($end > $#ids);
    } # while ($start <= $#ids)
    close($outfh) if ($output_isa_file);

    return $output;
}

# Bioperlize options
sub _hashref_to_bioperl {
    my $hash = shift || return;
    my %ret  = ();
    foreach my $key (keys %$hash) {
	next unless (defined $hash->{$key});
	next if (grep { $key eq $_ } @NON_BIOPERL_OPTIONS);
	my $nkey = ($key =~ /^-/) ? $key : "-$key";
	$ret{$nkey} = $hash->{$key};
    }
    return %ret;
}

1;
