# POD documentation - main docs before the code

=head1 NAME

Rotifer::DB::NCBI - interface to NCBI programs

=head1 SYNOPSIS

  # Creating a new parser
  use Rotifer::DB::NCBI;
  my @oldGIs = Rotifer::DB::NCBI->id2history($gi);

=head1 DESCRIPTION

Rotifer::DB::NCBI aggregates simplified interfaces for tools
developed at the NCBI.

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::DB::NCBI;

use strict;
use warnings;
use Carp::Clan qr/^Rotifer::DB::NCBI/;
use File::Basename;
use File::Find;
use File::Which qw(which);
use File::Monitored qw(:temp);
use List::MoreUtils qw(uniq);
use IO::File;
use IO::String;
use Rotifer::DB::NCBI::EUtils qw(get_seqids);
use Rotifer::DB::NCBI::EDirect;
use Rotifer::Utils qw(nr2ids find_file_by_name);
use Scalar::Util qw(reftype blessed);
use Sub::Exporter -setup => {
    exports => [qw/fastacmd get_replacements id2fasta id2history id2update id2dbxref protein2dna/],
};

=head2 get_replacements

 Title   : get_replacements
 Usage   : Rotifer::DB::NCBI->get_replacements(\%opts, @gis)
 Function: find/download nucleotide entries for protein identifiers
 Returns : list of IO::File objects
 Args    : array of sequence identifiers (strings)

=cut

sub get_replacements {
    my ($opts, @ids) = @_;

    # Esummary
    my %replaced = ();
    my %looks_like_acc = ();
    my @docsums = Rotifer::DB::NCBI::EUtils->esummary({ db => 'protein', debug => $opts->{debug}, retmax => $opts->{retmax} }, @ids);
    foreach my $docsum (@docsums) {
	my ($status) = $docsum->get_contents_by_name('Status');
	next if ($status eq 'live');
	my $id = $docsum->get_id;
	my ($nid) = $docsum->get_contents_by_name('ReplacedBy');
	if (defined $nid && length $nid) {
	    if ($nid =~ /^\d+$/) {
		$replaced{$nid} = $id;
	    } else {
		$nid =~ s/\.\d+$//;
		$looks_like_acc{$nid} = $id;
	    }
	}
    }

    # Accession to GI
    foreach my $arrayref (get_seqids(keys %looks_like_acc)) {
	my $origID = $looks_like_acc{$arrayref->[1]} || $looks_like_acc{$arrayref->[0]};
	next unless (defined $origID);
	$replaced{$arrayref->[0]} = $origID;
    }

    return %replaced;
}

=head2 protein2dna

 Title   : protein2dna
 Usage   : Rotifer::DB::NCBI->protein2dna(\%opts, @gis)
 Function: find/download nucleotide entries for protein identifiers
 Returns : list of IO::File objects
 Args    : array of sequence identifiers (strings)

 Users may precede the list of identifiers with an
 optional anonymous hash reference that may include
 any key supported by Rotifer::DB::NCBI::EUtils plus

 output => IO::File object
 debug => print progress messages for debugging
 retry => numbers of attempts to download the data

 gbksuffix => list of GenBank filename extensions
 dnasource => directories with nucleotide GenBank files
 proteinsource => directories with protein GenBank files

 !!!!!!WARNING!!!!!!WARNING!!!!!!WARNING!!!!!!

 Although this method tries to avoid downloading
 the same Genbank entry twice it does not yet cache
 requests to avoid repeated downloads between calls.

 Also, while caring to avoid multiple downloads, it
 is likely that the input order will be LOST and no
 effort is made to preserve it.

=cut

sub protein2dna {
    my $self = shift if (defined $_[0] && ($_[0] eq __PACKAGE__ || ref($_[0]) eq __PACKAGE__));
    my $opts = { %{ shift @_ } } if (defined $_[0] && UNIVERSAL::isa($_[0],"HASH"));

    # Process options
    $opts->{retry} = 5 unless (exists $opts->{retry});
    $opts->{debug} = 0 unless (exists $opts->{debug});

    # Accession to GI
    my %ids = ();
    map { $ids{$_} = 1 } grep { !/^\d+$/ } @_; # Accessions
    map { $ids{$_->[1]} = 1 } get_seqids(grep { /^\d+$/ } @_); # GI to accession.version
    #my %replaced = get_replacements($opts, keys %ids);
    return undef if (!scalar(keys %ids)); # No ids?

    # Output
    my $outfh = delete $opts->{output};
    if (!defined $outfh) {
	$outfh = IO::File->new_tmpfile;
    } elsif (!blessed $outfh) {
	if (reftype($outfh) eq 'GLOB') {
	    $outfh = IO::Handle->new_from_fd($outfh,"w");
	} else {
	    $outfh = IO::File->new($outfh,"w");
	}
    }
    #print STDERR join(" ","P2D",$outfh,@ids),"\n";

    # If GIs are missing, use elink to find related nucleotide sequences
    #my $missing = _fetch_nucleotide_using_ipg($opts, $outfh, \%ids);
    my $missing = _fetch_nucleotide_using_elink($opts, $outfh, \%ids); # ELink
    $missing = _fetch_nucleotide_using_efetch($opts, $outfh, \%ids) if ($missing); 
    #$missing = _fetch_nucleotide_using_idfetch($opts, $outfh, \%ids) if ($missing);

    # Rewind (if possible) and complain about things that didn't work
    $outfh->seek(0,0) if ($outfh->can("seek"));
    carp join(" ","protein2dna -> Server error, missing",$missing,":",sort keys %ids) if ($opts->{debug} && scalar(keys %ids));
    return ($outfh);
} # sub protein2dna

=head2 id2fasta

 Title   : id2fasta
 Usage   : Rotifer::DB::NCBI->id2fasta(\%opts, @gis)
 Function: use fastacmd and or id2fasta to retrieve sequences
 Returns : filename or file handle
 Args    : (optional) hash reference and 
           array of sequence identifiers (strings)

=cut

sub id2fasta {
    my $self = shift if (defined $_[0] && ($_[0] eq __PACKAGE__ || UNIVERSAL::isa($_[0],__PACKAGE__)));
    my $opts = { %{shift(@_)} } if (defined $_[0] && UNIVERSAL::isa($_[0],"HASH"));
    my @ids = @_;
    return "" unless (scalar @ids);

    my %gi = map { @{$_} } get_seqids(grep { /^\d+$/o } @ids);
    my %missing = ();
    for (my $i=0; $i<=$#ids; $i++) {
	$ids[$i] = $gi{$ids[$i]} if (exists $gi{$ids[$i]});
	$missing{$ids[$i]} = 1;
    }

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

    # Try EUtils
    if (scalar(keys %missing)) {
	my $fasta = ""; 
	my $io = IO::String->new($fasta);
	$io->setpos(length $fasta); # go to the end
	my $eopts = { %$opts, db => undef, output => $io };
	Rotifer::DB::NCBI::EUtils->id2fasta($eopts, keys %missing);
	foreach my $header (split(/\n/,$fasta)) {
	    next if ($header !~ /^>/);
	    foreach my $hash (nr2ids($header)) {
		delete $missing{ $hash->{accession} };
		delete $missing{ $hash->{accession} .".". $hash->{version} };
	    }
	}
	$io->close;
	print $outfh $fasta if (length $fasta);
    } # if (scalar(keys %missing))

    # Try id1_fetch if anything else fails and there are still some missing IDs: no full header though...
    if (scalar(keys %missing) && which("id1_fetch")) {
	carp join(" ","Retrieveing",scalar(@ids),"sequences using idfetch") if ($opts->{debug});

	# Temporary file
	my ($fh, $file) = ttempfile();
	print $fh join("\n",uniq sort keys %missing),"\n";
	close($fh);

	# Call idfetch
	carp "id1_fetch -fmt fasta -db protein -in $file" if ($opts->{debug});
	open(ID2FASTA, "id1_fetch -fmt fasta -db protein -in $file |") || croak "id1_fetch call failed!";
	my $load = 0;
	my $last = undef;
	while (<ID2FASTA>) {
	    if (/^>/) { # Is a header
		$load = 0; # Reset load controller
		my @hash = grep { exists $missing{$_->{accession}} || exists $missing{$_->{accession}.".".$_->{version}} } nr2ids($_);
		if (scalar @hash) {
		    map {
			delete $missing{$_->{accession}};
			delete $missing{$_->{accession} .'.'. $_->{version}};
		    } @hash;
		    $load = 1;
		}
		$_ = "\n$_" if ($opts->{clean} && defined $last && $last !~ /\n$/);
	    } else {
		s/\s//g if ($opts->{clean}); # Clean FASTA sequence
	    }
	    if ($load) {
		print $outfh $_;
		$last = $_;
	    }
	}
	close(ID2FASTA);

	print $outfh "\n" if (defined $last && $last =~ /\n$/);
    } #     if (scalar(keys %missing) && which("id1_fetch"))

    # Retrieve sequences from database using blastdbcmd
    if (scalar(keys %missing)) {
	carp join(" ","Retrieving",scalar(@ids),"sequences using blastdbcmd") if ($opts->{debug});
	my $exe = which("blastdbcmd");
	if (defined $exe) {
	    # Temporary file
	    @ids = map { s/^\s+//; s/\s+$//; $_ } @ids;
	    my ($fh, $file) = ttempfile();
	    print $fh join("\n",uniq sort keys %missing),"\n";
	    close($fh);

	    # blastdbcmd options
	    my @fopts = ();
	    if (exists $opts->{blastdbcmd_opts}) { # (key, value) to key only
		# Enable by default
		foreach my $flag (qw(-target_only -ctrl_a)) {
		    push(@fopts,$flag) if (!exists $opts->{blastdbcmd_opts}{$flag} || $opts->{blastdbcmd_opts}{$flag});
		    delete $opts->{blastdbcmd_opts}{$flag};
		}
	    }
	    push(@fopts, 
		 (defined $opts && exists $opts->{blastdbcmd_opts} ? %{ $opts->{blastdbcmd_opts} } : ()),
		 -db => $opts->{db} || "nr", 
		 -entry_batch => $file,
		);

	    carp "$exe @fopts 2> /dev/null\n" if ($opts->{debug});
	    open(BLASTDBCMD,"$exe @fopts 2> /dev/null |") || croak "Blastdbcmd call failed!";
	    my $load = 0;
	    my $last = undef;
	    while (<BLASTDBCMD>) {
		if (/^>/) {
		    $load = 0; # Reset load controller
		    my @hash = grep { exists $missing{$_->{accession}} || exists $missing{$_->{accession} .'.'. $_->{version}} } nr2ids($_);
		    if (scalar @hash) {
			map {
			    delete $missing{$_->{accession}};
			    delete $missing{$_->{accession} .'.'. $_->{version}};
			} @hash;
			$load = 1;
		    }
		    $_ = "\n$_" if ($opts->{clean} && defined $last && $last !~ /\n$/);
		} elsif ($opts->{clean}) {
		    s/\s//g; # Clean FASTA
		}
		if ($load) {
		    print $outfh $_;
		    $last = $_;
		}
	    }
	    close(BLASTDBCMD);
	    unlink($file);

	    print "\n" if (defined $last && $last !~ /\n$/);
	} # if (which("blastdbcmd"))
    } # if (scalar(keys %missing))

    close($outfh) if ($output_isa_file);
    return $output;
}

=head2 fastacmd

 Title   : fastacmd
 Usage   : Rotifer::DB::NCBI->fastacmd(\%opts, @gis)
 Function: use fastacmd to extract data from a BLAST database
 Returns : string (all sequences in FASTA format)
 Args    : (optional) hash reference and 
           array of sequence identifiers (strings)

=cut

sub fastacmd {
    my $self = shift if (defined $_[0] && ($_[0] eq __PACKAGE__ || UNIVERSAL::isa($_[0],__PACKAGE__)));
    my $opts = shift if (defined $_[0] && UNIVERSAL::isa($_[0],"HASH"));
    my %ids  = (); my @ids = map { $ids{$_} = 1; $_ } @_;
    return "" unless (scalar @ids);

    # Temporary file
    my ($fh, $file) = ttempfile();
    print $fh join("\n",sort keys %ids),"\n";
    close($fh);

    # Retrieve sequences from database
    if (exists $opts->{db} && defined $opts->{db}) { # Efetch compatibility
	$opts->{"-d"} = $opts->{db} eq "protein" ? "nr" : $opts->{db} =~ /nucleotide|nuccore/ ? "nt" : $opts->{db};
    } else {
	$opts->{"-d"} = "nr";
    }
    $opts->{"-i"} = $file;
    my @opts = %$opts;
    open(FASTACMD,"fastacmd @opts 2> /dev/null |");
    my $text = ""; 
    while (<FASTACMD>) {
	s/\n//g; # Clean FASTA
	$_ = "$_\n" if (/^>/);
	$_ = "\n$_" if (/^>/ && length $text);
	$text .= $_;
    } 
    close(FASTACMD);

    return $text;
}

=head2 id2update

 Title   : id2update
 Usage   : Rotifer::DB::NCBI->id2update(@gis)
 Function: get the newest GI using id1_fetch
 Returns : hash of integers (input GI is key)
 Args    : NCBI GIs or accessions

=cut

sub id2update {
    my $self = shift if (ref($_[0]) eq __PACKAGE__ || $_[0] eq __PACKAGE__);
    my %newest = id2history(@_);
    map { $newest{$_} = $newest{$_}->[0]{'accession'} } keys %newest;
    return %newest;
}

=head2 id2history

 Title   : id2history
 Usage   : Rotifer::DB::NCBI->id2history(@gis)
 Function: retrieve GI history using idfetch
 Returns : hash of arrays of hashes

           Returned hash structure: 
           KEYS   => input GIs/accessions (if found)
           VALUES => arrays of hashes with the keys

           history, accession, acctype, dbname, and version

           The rank values are for ordering: lower
           ranks are newer entries

 Args    : NCBI GIs or accessions

 Note    : any input not found in NCBI databases
           returns nothing (i.e. no key in the hash)

=cut

sub id2history {
    my $self = shift if (ref($_[0]) eq __PACKAGE__ || $_[0] eq __PACKAGE__);
    my @args = @_;
    @args = map { s/^\s+//; s/\s+$//; $_ } @args; # Erase trailing/leading spaces

    # Copy input to temporary file
    my ($fh,$file) = ttempfile();
    print $fh join("\n",@args),"\n";
    close($fh);

    # Call idfetch
    my $arrayref = undef; # Each history will be an anonymous array
    my %history = (); # This hash will include all GIs, old and new, as keys
    open(IDFETCH, "idfetch -i 4 -t 5 -G $file |") || croak "idfetch call failed!";
    while (<IDFETCH>) {
	if (/^GI/) { # New header: prepare to create a new history
	    $arrayref = undef;
	    next;
	}
	my ($gi) = /^(\d+)/;
	next if (!defined $gi || exists $history{$gi});
	if (defined $arrayref) {
	    push(@$arrayref, $gi);
	} else {
	    $arrayref = [ $gi ];
	}
	$history{$gi} = $arrayref;
    }
    close(IDFETCH);

    # Retrieve all ids
    my %out = ();
    my %dbxref = id2dbxref(@args, keys %history);
    foreach my $input (@args) {
	my $gi = $input;
	if (!exists $history{$input}) { # Accession?
	    if (exists $dbxref{$input}) {
		$gi = $dbxref{$input}->[0]{"accession"}; # use GI!
	    } else {
		next; # give up! :-(
	    }
	}
	foreach my $old (@{$history{$gi}}) {
	    foreach my $xref (@{$dbxref{$old}}) {
		push(@{$out{$input}}, { history => $old, %$xref });
	    }
	}
    }

    return %out;
}

=head2 id2dbxref

 Title   : id2dbxref
 Usage   : Rotifer::DB::NCBI->id2dbxref(@IDs)
 Function: find the GIs of any NCBI sequence ID
 Returns : hash of arrays of hashes

           The returned hash contains

           KEYS   => input IDs       (only if found in NCBI)
           VALUES => array of hashes (see Rotifer::Utils::nr2ids)

 Args    : array of NCBI GIs and accessions numbers

 Note    : queries that are not NCBI IDs will not be
           in the returned hash

=cut

sub id2dbxref {
    my $self = shift if (ref($_[0]) eq __PACKAGE__ || $_[0] eq __PACKAGE__);
    my @args = @_; @args = map { s/^\s+//; s/\s+$//; $_ } @args;
    my %args = map { ($_,undef) } @args;

    # Copy input to temporary file
    my ($fh,$file) = ttempfile();
    print $fh join("\n",@args),"\n";
    close($fh);

    # Call idfetch
    my %seen   = ();
    my @seqids = ();
    open(SEQIDS, "idfetch -i 2 -t 5 -G $file |") || croak "idfetch call failed!";
    while (<SEQIDS>) {
	chomp;
	next if (exists $seen{$_}); # Remove redundancy from equivalent accessions and GIs
	my $aoh = [ nr2ids($_) ];
	foreach my $dbxref (@$aoh) {
	    delete $dbxref->{accgroup};
	    delete $dbxref->{description};
	    $args{$dbxref->{"accession"}} = $aoh if (exists $args{$dbxref->{"accession"}});
	}
	$seen{$_} = 1;
    }
    close(SEQIDS);

    # Remove queries
    map { delete $args{$_} unless (defined $args{$_}) } @args;

    return %args;
}

# Downloading DNA sequences and checking for the presence of targets
#
# ELink returns related DNA entries that do not contain input protein GI
# thus forcing us to minimize network load by downloading one at a time
# and checking whether the GI is there
sub _fetch_nucleotide_using_elink {
    my ($opts, $outfh, $ids) = @_;

    # Search for corresponding nucleotide entries
    carp join(" ","processing",scalar(keys %$ids),"IDs using ELink/EFetch") if ($opts->{debug});
    my @map = Rotifer::DB::NCBI::EUtils->elink({ db       => 'nucleotide',
						 dbfrom   => 'protein',
						 debug    => $opts->{debug},
						 retmax   => $opts->{retmax},
					       }, keys %$ids);
    return scalar(keys %$ids) unless (scalar @map);

    # Load nucleotide docsums
    my %nsum = ();
    foreach my $docsum (Rotifer::DB::NCBI::EUtils->esummary({ db       => 'nucleotide',
							      debug    => $opts->{debug},
							      retmax   => $opts->{retmax},
							    }, map { $_->[1] } @map)) {
	$nsum{$docsum->get_id} = $docsum;
    }

    #  Sort nucleotide sequences by length
    @map = sort { 
	my ($alen) = $nsum{$a->[1]}->get_contents_by_name('Length') if (exists $nsum{$a->[1]});
	$alen = 0 if (!defined $alen);
	my ($blen) = $nsum{$b->[1]}->get_contents_by_name('Length') if (exists $nsum{$b->[1]});
	$blen = 0 if (!defined $blen);
	$a->[0] <=> $b->[0] || $blen <=> $alen || $a->[1] <=> $b->[1]
    } @map;

    # Check for matches to queries and print
    foreach my $map (@map) {
	next if (!exists $ids->{$map->[0]});
	next if ($map->[2] =~ /protein_nuccore_(mrna|wgs)/); # Avoid mRNA and WGS
	my $string = "";
	my $io = IO::String->new($string);
	my $fopts = { db => "nucleotide", output => $io, debug => $opts->{debug} };
	Rotifer::DB::NCBI::EUtils->id2genbank($fopts, $map->[1]);
	my @gi     = ($string =~ /\/db_xref=\"GI:(\d+)\"/g);
	my @accver = ($string =~ /\/protein_id=\"(\S+)\"/g);
	my $found = 0;
	foreach my $id (@gi,@accver) {
	    if (exists $ids->{$id}) { 
		delete $ids->{$id};
		$found = 1;
	    }
	}
	print $outfh $string if ($found);
    }

    return scalar(keys %$ids);
}

sub _fetch_nucleotide_using_efetch {
    my ($opts, $outfh, $ids) = @_;
    my @gbksuffixes = exists($opts->{gbksuffixes}) && scalar(@{$opts->{gbksuffixes}}) ? @{$opts->{gbksuffixes}} : qw(.gbk .gb .gbff .gbff.gz);
    require Rotifer::Parser;
    my $gb2hash = Rotifer::Parser->create("genbank2hashes", { annotation => [ 'all' ], include => [ qw/CDS/ ] });

    carp join(" ","processing",scalar(keys %$ids),"IDs using EFetch protein/parse/EFetch nucleotide/drop RNA") if ($opts->{debug});
    my $fopts = { rettype => 'gbwithparts', retmode => 'text', debug => $opts->{debug} };

    # Parse the protein's GenBank file 
    my %dnaids  = ();
    my %dnafiles = map { my @g = fileparse($_,@gbksuffixes); ($g[0],$_) } @{$opts->{dnasource}};
    foreach my $pfile (@{$opts->{proteinsource}},"Rotifer::DB::NCBI::EUtils") {
	my $protfile = $pfile;

	# Download if needed
	my $istemp = 0;
	if ($protfile eq "Rotifer::DB::NCBI::EUtils") {
	    carp join(" ","processing",scalar(keys %$ids),"IDs using EFetch protein") if ($opts->{debug} > 5);
	    $protfile = Rotifer::DB::NCBI::EUtils->id2genbank($fopts, keys %$ids);
	    $istemp = 1;
	}

	# Parse
	my ($header, $keys, @table) = $gb2hash->parse($protfile);
	foreach my $hash (@table) {
	    next unless (exists $hash->{coded_by});
	    my $acc = $hash->{coded_by};
	    $acc =~ s/complement\((.+)\)/$1/;
	    $acc =~ s/:\S+//;
	    next if ($acc =~ /^(N|X)M_\d+/);
	    if ($acc !~ /^[A-Z]/) {
		carp "Strange accession $acc";
		next;
	    }
	    $dnaids{$acc}++ unless (exists $dnafiles{$acc});
	}

	unlink($protfile) if ($istemp);
    } # foreach my $protfile (@protfile)
    %dnaids = scalar(%dnaids) ? map { ($_->[0],$_->[1]) } get_seqids(keys %dnaids) : ();

    # Downloading GenBank file for DNA sequences using accessions
    foreach my $gbfile (@{$opts->{dnasource}},"Rotifer::DB::NCBI::EUtils") {
	my $dnafile = $gbfile;

	# Open DNA file
	my $istemp = 0;
	my $dnafh = undef;
	if ($dnafile =~ /\.gz$/) {
	    open($dnafh,"gunzip -c $dnafile |");
	} else {
	    if ($dnafile eq "Rotifer::DB::NCBI::EUtils") {
		($dnafile) = Rotifer::DB::NCBI::EUtils->id2genbank($fopts,sort { $a <=> $b } keys %dnaids);
		$istemp = 1;
	    }
	    open($dnafh,"<$dnafile");
	}

	# Parse DNA file
	my $in = 0;
	while (<$dnafh>) {
	    /LOCUS/ && do {
		$in = /\b(mRNA|RNA)\b/ ? 0 : 1;
	    };

	    /VERSION +(\S+)/ && do {
		delete $dnaids{$1};
		/GI:(\d+)/ && do { delete $dnaids{$2} };
	    };

	    if ($in) {
		if (/^\s*\/db_xref=\"GI:(\d+)\"/) {
		    delete $ids->{$1};
		} elsif (/\s*\/protein_id=\"(\S+)\"/) {
		    delete $ids->{$1};
		}
		print $outfh $_;
	    }

	    /^\/\// && do {
		$in = 0;
	    };
	}
	close($dnafh);

	unlink($dnafile) if ($istemp);
    } # if (scalar(keys %dnaids))

    return scalar(keys %$ids); # Updating list
}

sub _fetch_nucleotide_using_idfetch {
    my ($opts, $outfh, $ids) = @_;

    # Choose NCBI tool program
    my ($method,$program) = ();
    if (($program) = which("edirect.pl")) {
	$method = "$program -fetch -format gbwithparts -db nucleotide -id";
    } elsif (($program) = which("idfetch")) {
	$method = "$program -t 3 -g";
    } elsif (($program) = which("id1_fetch")) {
	$method = "$program -fmt genbank -gi";
    } else {
	carp "NCBI tools (id1_fetch or idfetch) for downloading DNA data are not available in your system!!!!
Please install NCBI's C or C++ toolkit!!!
Falling back to the much slower EUtils-only approach...";
    }

    # Use $method to download DNA Genbank flat files directly
    if (defined $method && scalar(keys %$ids)) {
	carp join(" ","processing",scalar(keys %$ids),"IDs using $method: ",keys %$ids) if ($opts->{debug});
	foreach my $id (keys %$ids) { # Process one protein at a time to avoid repeated downloads
	    my $in = 0;
	    my $string = ""; # Load entire file into memory
	    my @found  = ();
	    open(my $id1, "$method $id |");
	    while (<$id1>) {
		next if (/^\s*$/);
		last if (/^LOCUS\s+.+\b(mRNA|RNA)\b/); # Ignore (m)RNA entries
		if (/protein_id=\"(\S+)\"/) {
		    if (exists $ids->{$1}) {
			push(@found,$1);
		    }
		}
		if (/db_xref=\"GI:(\d+)\"/) {
		    if (exists $ids->{$1}) {
			push(@found,$1);
		    }
		}
		$string .= $_;
	    } # while (<$id1>)
	    close($id1);
	    if (scalar @found && $string =~ /^LOCUS.+\/\/\s*$/s) { # Complete!
		print $outfh $string;
		map { delete $ids->{$_} } @found;
	    }
	} # foreach my $id (@ids)
    } # if (defined $method)

    return scalar(keys %$ids);
}

sub _fetch_nucleotide_using_ipg {
    my ($opts, $outfh, $ids) = @_;

    # Options
    my $debug  = $opts->{debug}  || 0;
    my $retmax = $opts->{retmax} || 2000;
    my $db     = $opts->{db}     || "protein";

    # Retrieve data for each sequence
    foreach my $id (keys %$ids) {
	# Retrieve and parse IPG report using EDirect's efetch
	my $ipg;
	efetch({ output => \$ipg }, $id);
	my @ipgs = split(/\n/,$ipg);

	# Parse Identical Protein Report
	my %header = ();
	my @dna = ();
#	while (<IPG>) {
#	    chomp;
#	    my @row = split(/\t/);

	    # Process header
#	    if (!scalar(keys %header)) {
#		%header = map { $_ } @row;
#		next;
#	    }

	    # Process data row
	    
#	}
    }
}

# Make perl compiler happy
1;
