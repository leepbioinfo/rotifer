# POD documentation - main docs before the code

=head1 NAME

Rotifer::DB::NCBI::EDirect - an API to access NCBI's EUtils services
                             using EDirect Utilities 

=head1 SYNOPSIS

  # Convert protein identifiers to DNA GI's
  use Rotifer::DB::NCBI::EDirect qw(elink);
  my @table = elink({ from => "protein", to => "nucleotide" }, $gi);

=head1 DESCRIPTION

Rotifer::DB::NCBI::EDirect provides a simple interfaces to retrieve
data from NCBI E-Utils service. It is built on top of NCBI's Edirect
utilities, a collection of pipe-friendly utilities maintained by NCBI
itself and integrated using IPC::Run.

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 DEPENDENCIES

=over

=item IPC::Run

=item Moo

=item EDirect utilities

=back

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::DB::NCBI::EDirect;

use strict;
use warnings;
use autodie qw(:all);
use Carp::Clan qr/^Rotifer::DB::NCBI::EDirect/;
use File::Monitored qw(ttempfile);
use IO::String;
use IO::File;
use IPC::Run qw(run start pump finish timeout);
use Scalar::Util qw(blessed reftype);
use Sub::Exporter -setup => {
    exports => [qw/id2genbank id2fasta
                   get_seqids accession2gi gi2accession
                   elink efetch epost esummary/],
};

# Global variables
our $MAX_DOWNLOAD_ATTEMPTS = 3;
our $DEFAULT_RETMAX        = 2000000;
#our $DEFAULT_SERVER        = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils';
#our %DEFAULT_SERVICES      = ( 'elink' => 'elink.fcgi' );

=head2 id2fasta

 Title   : id2fasta
 Usage   : Rotifer::DB::NCBI::EDirect->id2fasta(\%opts, @gis)
 Function: retrieve sequences in FASTA format
 Returns : file name or file handle (same as key "output")
 Args    : array of sequence identifiers (strings)

 Users may precede the list of identifiers with an
 optional anonymous hash reference that may include
 any option suuported by 

   debug => print progress messages for debugging
  output => file name or file handle
   retry => numbers of attempts to download the data

=cut

sub id2fasta {
    my $self = shift if (defined $_[0] && ($_[0] eq __PACKAGE__ || ref($_[0]) eq __PACKAGE__));
    my %opts = %{shift(@_)} if (defined $_[0] && ref($_[0]) eq "HASH");
    my @ids  = @_;

    # Process options
    my $attempts = $opts{retry} || $MAX_DOWNLOAD_ATTEMPTS;
    my $debug    = $opts{debug} || 0;
    $opts{'-format'} = 'fasta';
    $opts{'-mode'}   = 'text';

    # Ouput
    my ($outfh,$output);
    my $output_isa_file = 1;
    if (exists $opts{output} && defined $opts{output}) {
	$output = delete $opts{output};
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

    # Query, retrieve, process...
    my $round  = 0;
    my %ids = map { ($_,1) } @ids;
    my $last = "";
    while (scalar(@ids) && $round < $attempts) {
	carp join(" ","id2fasta -> Server error, missing",scalar(@ids),". Redo #$round:",@ids) if ($debug > 1 && $round > 0);
	carp join(" ","id2fasta -> Retrieving",scalar(@ids),"of",scalar(@_),"in round #",$round) if ($debug);
	my $filename = efetch(\%opts, @ids);
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
		s/\s//g if ($opts{clean});
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
 Usage   : Rotifer::DB::NCBI::EDirect->id2genbank(\%opts, @gis)
 Function: download genbank files
 Returns : file name or undef
 Args    : array of sequence identifiers (strings)

 Users may precede the list of identifiers with an
 optional anonymous hash reference that may include
 options supported by Rotifer::DB::NCBI::EDirect::efetch

   debug => print progress messages for debugging
  output => file name or file handle to dump output
   retry => numbers of attempts to download the data

=cut

sub id2genbank {
    my $self = shift if (defined $_[0] && ($_[0] eq __PACKAGE__ || ref($_[0]) eq __PACKAGE__));
    my %opts = %{shift(@_)} if (defined $_[0] && ref($_[0]) eq "HASH");
    my @ids  = @_;

    # Process options
    my $retry = $opts{retry} || $MAX_DOWNLOAD_ATTEMPTS;
    my $debug = $opts{debug} || 0;
    $opts{"-format"} = 'gbwithparts';
    $opts{"-mode"}   = 'text';

    # Output
    my ($outfh, $output) = exists $opts{output} ? (undef,delete $opts{output}) : &ttempfile();
    if (!ref $output) {
	open($outfh,">$output") || die "Could not create outut file $output";
    } else {
	$outfh  = $output;
	$output = undef;
    }

    # GI to accession
    my %gi = map { @$_ } get_seqids(grep { /^\d+$/ } @_);
    if (scalar(keys %gi)) {
	for (my $i=0; $i<=$#ids; $i++) {
	    $ids[$i] = $gi{$ids[$i]} if (exists $gi{$ids[$i]});
	}
    }

    # Query, retrieve, process...
    my $round  = 0;
    while (scalar(@ids) && $round < $retry) {
	carp join(" ","id2genbank -> Server error, missing",scalar(@ids),". Redo #$round:",@ids) if ($debug > 1 && $round > 0);
	carp join(" ","id2genbank -> Retrieving",scalar(@ids),"of",scalar(@_),"in round #",$round) if ($debug);
	$opts{output} = IO::File->new_tmpfile; # Second temporary file: raw download
	efetch(\%opts, @ids);

	#++$round && next unless (defined $temp && -s $temp);
	++$round && next unless (defined $opts{output});

	# Fetching
	my %found  = ();
	my $string = "";
	#open (FETCHED,"<$temp") || die "Could not open temporary file $temp";
	#while (<FETCHED>) {
	while (<$opts{output}>) {
	    next if /^\s*$/;

	    # End of sequence: store
	    if (/^\/\/\s*$/) {
		if (length $string && $string =~ /^LOCUS/) {
		    $string .= $_;
		    print $outfh $string;
		}
		$string = "";
		%found = ();
		next;
	    }

	    # Accession
	    if (/^ACCESSION\s+(\S+)/) {
		$found{$1} = 1;
	    }

	    # Accession.Version
	    elsif (/^VERSION\s+(\S+)/) {
		$found{$1} = 1;
		if (/GI:(\d+)/) { # GI
		    $found{$1} = 1;
		}
	    }

	    $string .= $_; # Concatenate Genbank file
	} # while (<$opts{output}>)
	close($opts{output});
	delete $opts{output};

	@ids = grep { exists $found{$_} } @ids;
	$round++;
    } # while  (scalar(@ids) && $round < $retry)
    close($outfh) if (defined $outfh);

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
 Args    : array of XML strings returned by epost

=cut

sub get_count {
    my $self = shift if (defined $_[0] && ($_[0] eq __PACKAGE__ || ref($_[0]) eq __PACKAGE__));

    my $sum = 0;
    foreach my $history (@_) {
	my @count = ($history =~ /<Count>(\d+)<\/Count>/);
	map { $sum += $_ } @count;
    }

    return $sum;
}

=head2 get_seqids

 Title   : get_seqids
 Usage   : Rotifer::DB::NCBI::EDirect->get_seqids(@ids)
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

    # Options
    my $db     = $opts->{db}     || "protein";
    my $debug  = $opts->{debug}  || 0;
    my $retmax = $opts->{retmax} || 200;
    my $retry  = $opts->{retry}  || $MAX_DOWNLOAD_ATTEMPTS;

    # Send/retrieve
    my @array = ();
    my $round = 0;
    while (scalar(@ids) && $round < $retry) { # Shit happens... try again!
	carp join(" ","missing",scalar(@ids),". Redo #$round:",@ids) if ($debug > 1 && $round > 0);
	carp join(" ","retrieving ",scalar(@ids),"of",scalar(@_),"in round #",$round) if ($debug);

	# Fetching
	my %found = ();
	my $start = 0;
	my $end   = $#ids > $retmax-1 ? $retmax-1 : $#ids;
	while ($start <= $#ids) {
	    # Call efetch
	    carp join(" ","sending chunk ${start}..$end of",scalar(@ids)) if ($debug);
	    my ($out, $err);
	    my $in = join(",",@ids[$start..$end]);
	    run [ "efetch", "-format", "seqid", "-db", $db, "-id", $in  ], \undef, ">pipe", \*FETCHED || carp $?;

	    # Parse seqids
	    my $acc = undef;
	    my $version = undef;
	    while (<FETCHED>) { # Parse ASN.1
		next if /^\s*$/;
		if (/accession\s+\"(\S+)\"/) {
		    $acc = $1;
		    $found{$acc} = 1;
		} elsif (/version\s+(\d+)/) {
		    $version = $1;
		} elsif (/Seq-id\s+::=\s+gi\s+(\d+)/) {
		    next unless (defined $acc);
		    my $gi = $1;
		    $found{$gi} = 1;
		    $acc .= ".$version" if (defined $version);
		    $found{$acc} = 1;
		    push(@array, [ $gi, $acc ]);
		    undef($version);
		    undef($acc);
		};
	    }

	    # Prepare for next batch
	    $start = $end+1;
	    $end   = $end+$retmax;
	    $end   = $#ids if ($#ids < $end);
	    close(FETCHED);
	} # while ($#ids > $start)

	@ids = grep { !exists $found{$_} } @ids;
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
 Usage   : Rotifer::DB::NCBI::EDirect->esummary(\%opts, @gis)
 Function: use the esummary service
 Returns : array of XML strings
 Args    : array of sequence identifiers (strings)

 Users may precede the list of identifiers with an
 optional anonymous hash reference that may include
 options supported by Rotifer::DB::NCBI::EDirect::efetch

=cut

sub esummary {
    my $self = shift if (defined $_[0] && ($_[0] eq __PACKAGE__ || ref($_[0]) eq __PACKAGE__));
    my %opts = defined ($_[0]) && ref($_[0]) eq "HASH" ? %{ shift(@_) } : {};

    # Options
    my $output;
    my $debug = $opts{debug} || 0;
    $opts{"-format"} = "docsum";
    $opts{"-mode"}   = "xml";
    $opts{"output"}  = \$output;
    efetch(\%opts,@_);
    return $output;

    # Epost: send ID list
    my @history = epost(\%opts, @_) if (scalar @_);
    return () if (!scalar @history);

    # Retrieve summaries
    my @docsum = ();
    foreach my $history (@history) {
	carp join(" ","sending next chunk of",scalar(@_)) if ($debug);
	my ($out, $err);
	run([ "esummary" ], \$history, \$out, \$err);
	carp $err if ($debug && defined $err && length $err);
	push(@docsum, $out) if  (defined $out && length $out);
    }

    return @docsum;
}

=head2 epost

 Title   : epost
 Usage   : Rotifer::DB::NCBI::EDirect->epost(\%opts, @gis)
 Function: use the epost service
 Returns : array of xml text containing history identifiers
 Args    : array of sequence identifiers (strings)

 Users may precede the list of identifiers with an
 optional anonymous hash reference that may include
 any of the following keys:

       db => database name
    debug => print progress messages for debugging
   retmax => maximum number of queries per download

=cut

sub epost {
    my $self = shift if (defined $_[0] && ($_[0] eq __PACKAGE__ || ref($_[0]) eq __PACKAGE__));
    my $opts = shift if (defined $_[0] && ref($_[0]) eq "HASH");
    my @ids = @_;

    # Options
    my $debug  = $opts->{debug}  || 0;
    my $retmax = $opts->{retmax} || $DEFAULT_RETMAX;
    my $db     = $opts->{db}     || "protein";

    # GI to accession
    my %gi = map { @$_ } get_seqids(grep { /^\d+$/ } @_);
    if (scalar(keys %gi)) {
	for (my $i=0; $i<=$#ids; $i++) {
	    $ids[$i] = $gi{$ids[$i]} if (exists $gi{$ids[$i]});
	}
    }

    # Epost: send list
    my @xml     = ();
    my $start   = 0;
    my $end     = $#ids > $retmax-1 ? $retmax-1 : $#ids;
    while ($start <= $#ids) {
	carp join(" ","Sending chunk ${start}..$end of",scalar(@ids),":",@ids) if ($debug);
	my ($out, $err);
	my $in = join("\n",@ids[$start..$end]);
	run([ "epost", "-db", $db, "-format", "acc" ], \$in, \$out, \$err);
	carp $err if ($debug && defined $err && length $err);
	push(@xml, $out) if  (defined $out && length $out);
	$start = $end+1;
	$end   = $end+$retmax;
	$end   = $#ids if ($#ids < $end);
    }

    return @xml;
}

=head2 elink

 Title   : elink
 Usage   : Rotifer::DB::NCBI::EDirect->elink(\%opts, @gis)
 Function: use the elink service
 Returns : array of XML text with history identifiers
 Args    : array of sequence identifiers (strings)

 Users may precede the list of identifiers with an
 optional anonymous hash reference with keys

    debug => print progress messages for debugging

 Other options must be supported by EDirect's elink
 and, therefore, start with a dash ("-").
 Example:

 elink({ -name => "protein_structure" }, @ARGV)

 Default values for required ELink options are:

    db     => protein
    target => nucleotide

 See

 http://eutils.ncbi.nlm.nih.gov/entrez/query/static/entrezlinks.html

 for a list of link names

=cut

sub elink {
    my $self = shift if (defined $_[0] && ($_[0] eq __PACKAGE__ || ref($_[0]) eq __PACKAGE__));
    my $opts = shift if (defined $_[0] && ref($_[0]) eq "HASH");

    # Options
    my $debug  = $opts->{debug}  || 0;
    my $db     = $opts->{dbfrom} || "protein";
    my $target = $opts->{db}     || "nucleotide";

    # Separate history and identifiers
    my @history = grep { /<ENTREZ_DIRECT>/ } @_; # History entries first!
    push(@history, epost($opts, grep { !/<ENTREZ_DIRECT>/ } @_));
    return () if (!scalar @history);

    # Calling Edirect's elink
    my @xml = ();
    foreach my $history (@history) {
	my ($out, $err);
	run([ "elink", "-target" ], \$history, \$out, \$err);
	carp $err if ($debug && defined $err && length $err);
	push(@xml, $out) if  (defined $out && length $out);
    }

    return @xml;
}

=head2 efetch

 Title   : efetch
 Usage   : Rotifer::DB::NCBI::EDirect->efetch(\%opts, @gis)
 Function: use the efetch service
 Returns : filehandle to iterate over results
 Args    : array of sequence identifiers (strings)
           or array of NCBI's XML history strings

 Users may precede the list of identifiers with an
 optional anonymous hash reference that may include:
 
 1. Any option supported by EDirect's efetch, provided
    its leading dash ("-") is preserved

 2. Any of the arguments below:

   debug => print progress messages for debugging
  output => file name or file handle to dump output
   epost => (boolean) whether to use epost (default 1)
   retry => numbers of attempts to download the data

 Note: when the optional hash's key 'output' is set to
       to a glob this subroutine will assume the value
       points to an open filehandle and will attempt to
       print all downloaded content to this filehandle.
       It will also return 'undef' instead of a file path.

=cut

sub efetch {
    my $self = shift if (defined $_[0] && ($_[0] eq __PACKAGE__ || ref($_[0]) eq __PACKAGE__));
    my $opts = shift if (defined $_[0] && ref($_[0]) eq "HASH");
    if (!exists $opts->{epost}) {
	if (scalar(@_) < 250) {
	    return _efetch_only($opts, @_);
	} else {
	    return _epost_efetch($opts, @_);
	}
    } elsif ($opts->{epost}) {
	return _epost_efetch($opts, @_);
    } else {
	return _efetch_only($opts, @_);
    }
}

=head1 Internal methods

=head2 _epost_efetch

 Title   : _epost_efetch
 Usage   : Rotifer::DB::NCBI::EDirect->_epost_efetch(\%opts, @gis)
 Function: epost then download with efetch
 Returns : IO::File object
 Args    : array of sequence identifiers (strings)
           or NCBI's XML history strings

 Users may precede the list of identifiers with an
 optional anonymous hash reference that may include:
 
 1. Any option supported by EDirect's efetch, provided
    its leading dash ("-") is preserved

 2. Any of the arguments below:

      db => database name, ignored if a key '-db' is found
   debug => print progress messages for debugging
  output => file name or file handle to dump output

=cut

sub _epost_efetch {
    my $self = shift if (defined $_[0] && ($_[0] eq __PACKAGE__ || ref($_[0]) eq __PACKAGE__));
    my %opts = %{ shift(@_) } if (defined $_[0] && ref($_[0]) eq "HASH");

    # Parse options that are not supported by EDirect 
    my $debug  = delete $opts{debug}  || 0;
    my $retmax = delete $opts{retmax} || $DEFAULT_RETMAX;
    my $output = delete $opts{output} || IO::File->new_tmpfile;
    $output = \$output unless (ref $output);
    $opts{"-db"} = (exists $opts{db} ? delete $opts{db} : "protein") unless (exists $opts{"-db"});
    $opts{"-format"} = "native" unless (exists $opts{"-format"});
    map { delete $opts{$_} } grep { !/^-/ } keys %opts;

    # Isolate histories, replace accessions with GIs and drop unknown accessions
    my (@history,@ids) = ();
    map { /<ENTREZ_DIRECT>/ ? push(@history,$_) : push(@ids,$_) } @_;
    #my $count = get_count(@history);

    # Epost: send GI list
    if (scalar @ids) {
	push(@history, epost({ db => $opts{"-db"}, debug => $debug, retmax => $retmax }, @ids));
	#$count += scalar(@ids);
    }
    return () if (!scalar @history);

    # Iterate over histories
    for (my $i=0; $i<=$#history; $i++) {
	carp join(" ","_epost_efetch -> Retrieving chunk",$i,"for",scalar(@ids),"queries.") if ($debug);
	run([ "efetch", %opts ], \$history[$i], $output, "2>pipe", \*ERR);
	my $err = join("",<ERR>);
	if (length $err) {
	    carp "EDirect's efetch call produced an error:\n",$err;
	}
    }

    return $output;
}

=head2 _efetch_only

 Title   : _efetch_only
 Usage   : Rotifer::DB::NCBI::EDirect->_efetch_only(\%opts, @gis)
 Function: download with efetch
 Returns : path to file with all data or undef (see note below)
 Args    : array of sequence identifiers (strings)
           or NCBI's XML history strings

 Users may precede the list of identifiers with an
 optional anonymous hash reference that may include:
 
 1. Any option supported by EDirect's efetch, provided
    its leading dash ("-") is preserved

 2. Any of the arguments below:

      db => database name, ignored if a key '-db' is found
   debug => print progress messages for debugging
  output => file name or file handle to dump output

=cut

sub _efetch_only {
    my $self = shift if (defined $_[0] && ($_[0] eq __PACKAGE__ || ref($_[0]) eq __PACKAGE__));
    my %opts = %{ shift(@_) } if (defined $_[0] && ref($_[0]) eq "HASH");

    # Isolate histories
    my (@history, @ids);
    map { /<ENTREZ_DIRECT>/ ? push(@history,$_) : push(@ids,$_) } @_;

    # Parse options
    my $debug  = delete $opts{debug}  || 0;
    my $retmax = delete $opts{retmax} || 200;
    my $output = delete $opts{output} || IO::File->new_tmpfile;
    $output = \$output unless (ref $output);
    $opts{"-db"} = (exists $opts{db} ? delete $opts{db} : "protein") unless (exists $opts{"-db"});
    $opts{"-format"} = "native"  unless (exists $opts{"-format"});
    map { delete $opts{$_} } grep { !/^-/ } keys %opts;

    # Fetching
    my $start = 0;
    my $end   = $#ids > $retmax-1 ? $retmax-1 : $#ids;
    while ($start <= $#ids) {
	# Call efetch
	carp join(" ","_efetch_only -> retrieving chunk ${start}..$end of",scalar(@ids)) if ($debug);
	my ($err);
	my $in = join(",",@ids[$start..$end]);
	run([ "efetch", , "-id", $in, %opts  ], \undef, $output, \$err) || carp "ERROR _efetch_only:\n",$err;

	# Prepare for next batch
	$start = $end+1;
	$end   = $end+$retmax;
	$end   = $#ids if ($#ids < $end);
    } # while ($#ids > $start)

    return $output;
}

1;
