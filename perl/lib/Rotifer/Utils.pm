# POD documentation - main docs before the code

=head1 NAME

Rotifer::Utils - some helper methods for Rotifer

=head1 SYNOPSIS

  # Documenting plugins
  use Rotifer::Utils;
  print Rotifer::Utils->describe_subclasses("Rotifer::DBIC::IO"),"\n";

=head1 DESCRIPTION

Rotifer::Utils aggregates useful routines for various tasks.

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 DEPENDENCIES

=over

=item Module::Find

=back

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::Utils;

#use namespace::autoclean; <-- good for Moose, kills Sub::Exporter
use autodie qw(:all);
use Carp::Clan qr/^Rotifer::Utils/;
use Cwd;
use File::Basename;
use File::Find;
use List::MoreUtils qw(any);
use Scalar::Util qw(looks_like_number);
use Sub::Exporter -setup => {
    exports => [qw/describe_subclasses find_file_by_name ids2nr nr2ids
                   pod2summary ref2text ruler relative2absolute
                   array_of_arrays_to_text aoa2txt array_of_arrays_to_tsv aoa2tsv/],
};
use strict;
use warnings;

# CLASS VARIABLES

our %Sources = (
    dbj => "DDBJ",
    emb => "EMBL",
    gi  => "GI",
    gb  => "Genbank",
    gnl => "General",
    lcl => "LocalDB",
    pdb => "PDB",
    pir => "PIR",
    prf => "PRF",
    ref => "RefSeq",
    sp  => "SwissProt",
    tpe => "TPE",
    tpd => "TPD",
    tpg => "TPG",
    );
map { $Sources{uc $_} = $Sources{$_} } keys %Sources;
our %ReversedSources   = map { ($Sources{$_},$_) } keys %Sources;
our @SequenceDatabases = sort keys %ReversedSources;

our %SeqIdTypeRank = (
    GI            => 0,
    ACCESSION     => 1,
    fakeGI        => 2,
    LOCUS         => 3,
    UnknownID     => 4,
    GI_fragment   => 5
    );
our @SeqIdTypes = keys %SeqIdTypeRank;

=head2 describe_subclasses

 Title   : describe_subclasses
 Usage   : $full_path = describe_subclasses("Rotifer::DBIC::IO")
 Function: load, parse NAME sections for subclasses and print
 Returns : 
 Args    : upper class name

=cut

sub describe_subclasses {
    my $opts = shift if (defined $_[0] && UNIVERSAL::isa($_[0],"HASH"));
    return unless (scalar @_);
    unshift(@_, $opts) if (defined $opts);
    my @arrayOfArrays = pod2summary(@_);
    my $header = [ 'Name', 'Description' ];
    push(@$header,grep { $_ ne "NAME" } @{$opts->{section}}) 
	if (defined $opts && exists $opts->{section});
    push(@$header, "Module");
    return array_of_arrays_to_text($header, @arrayOfArrays);
}

=head2 find_file_by_name

 Title   : find_file_by_name
 Usage   : find_file_by_name([ @dirs ], [ @names ], [ @suffixes ])
 Function: scan directories and return a list of matching files
 Returns : array of file paths
 Args    : three array references
           1 - list of directories
           2 - list of base file names
           3 - (optional) list of suffixes that must be matched

=cut

sub find_file_by_name {
    my $self = shift if (defined $_[0] && ($_[0] eq __PACKAGE__ || ref($_[0]) eq __PACKAGE__));
    my $opts = shift if (defined $_[0] && ref($_[0]) eq 'HASH');
    my ($directories,$basenames,$suffixes) = @_;
    my %bnames = map { ($_,1) } @$basenames;
    my @path = ();
    find({
	wanted => sub {
	    my ($target) = fileparse($_, @$suffixes);
	    push(@path, $File::Find::name) if (exists $bnames{$target} || exists $bnames{$_});
	}}, @$directories);
    return @path;
}

=head2 array_of_arrays_to_text

 Title   : array_of_arrays_to_text or aoa2txt
 Usage   : $txt = array_of_arrays_to_text([0,1], [1,0])
 Function: formats a matrix (array of arrays) for printing
 Returns : string
 Args    : array of arrays

 Aliases : aoa2txt

=cut

sub array_of_arrays_to_text {
    my $opts = shift if (defined $_[0] && ref($_[0]) eq 'HASH');

    # Calculate lengths
    my @length = ();
    foreach my $table (@_) {
	foreach my $i (0..$#{$table}) {
	    $length[$i] = 0 unless (defined $length[$i]);
	    my $str = defined $table->[$i] ? $table->[$i] : "";
	    $length[$i] = length($str) > $length[$i] ? length($str) : $length[$i];
	}
    }

    # Print
    my $line = "+-".join("-+-",map { "-" x $length[$_] } 0..$#length)."-+\n";
    my $text = "";
    my $j = 0;
    my $empty = exists $opts->{empty} ? $opts->{empty} : "";
    foreach my $row (@_) {
	$text .= "| ".join(" | ",map { sprintf("%-$length[$_]s",defined $row->[$_] ? $row->[$_] : $empty) } 0..$#length)." |\n";
	$text .= "$line" if ($j == 0 && $#_ > 0);
	$j++;
    }
    $text = "$line$text$line";

    return $text;
}

{no warnings 'once';
 *aoa2txt = \&array_of_arrays_to_text;
}

=head2 array_of_arrays_to_tsv

 Title   : array_of_arrays_to_tsv or aoa2tsv
 Usage   : $txt = array_of_arrays_to_tsv([0,1], [1,0])
 Function: output a matrix (array of arrays) as rows of
           text separated by TABS
 Returns : string
 Args    : array of arrays

 Aliases : aoa2tsv

=cut

sub array_of_arrays_to_tsv {
    my $opts = shift if (defined $_[0] && ref($_[0]) eq 'HASH');
    my $sep = $opts->{output_delimiter} || "\t";
    $opts->{empty} = "" unless (exists $opts->{empty});

    my $nofCols = -1;
    my @collen  = ();
    foreach my $row (@_) {
	$nofCols = $#{$row} if ($nofCols < $#{$row});
	if ($opts->{'pad'}) {
	    foreach my $col (0..$#{$row}) {
		if (!defined $collen[$col] || $collen[$col] < length($row->[$col] || "")) {
		    $collen[$col] = length($row->[$col] || "");
		}
	    }
	}
    }

    # Print
    my $text = "";
    foreach my $rowRef (@_) {
	my @row = ();
	for (my $i=0; $i<=$nofCols; $i++) {
	    my $elm = defined $rowRef->[$i] ? $rowRef->[$i] : $opts->{empty};
	    if ($opts->{'pad'}) {
		my $len = $collen[$i];
		if (exists $opts->{align} && $opts->{align}->{$i} eq 'right') {
		    $elm = sprintf("%${len}s",$elm);
		} else {
		    $elm = sprintf("%-${len}s", $elm);
		}
	    }
	    push(@row, $elm);
	}
	$text .= join($sep, @row)."\n";
    }

    return $text;
}

{no warnings 'once';
 *aoa2tsv = \&array_of_arrays_to_tsv;
}

=head2 pod2summary

 Title   : pod2summary
 Usage   : $full_path = pod2summary("Rotifer::DBIC::IO")
 Function: parse NAME sections from the POD documention
           of a file, class and subclasses
 Returns : array of arrays (three columns: basename, subclass, NAME) 
 Args    : (string) upper class name

=cut

sub pod2summary {
    my $opts = shift if (defined $_[0] && UNIVERSAL::isa($_[0],"HASH"));
    my @sections = qw(NAME);
    push(@sections, grep { $_ ne "NAME" } @{$opts->{section}}) 
	if (exists $opts->{section});

    my @desc = ();
    foreach my $class (@_) {
	# Load all classes from a path except the root class 
	use Module::Find;
	my @plugins = ($class, findallmod("$class")); # Module::Find in action

	# Parse each plugin
	foreach my $plugin (sort @plugins) {
	    eval "require $plugin" unless ( -f $plugin );
	    (my $name = $plugin) =~ s/.+:://;
	    (my $path = $plugin . ".pm") =~ s/::/\//g;

	    if ($@) {
		push(@desc, [ basename($name), " ===> won't compile??? I can't get a summary...", $plugin ]);
		next;
	    }

	    my @cols = ();
	    foreach my $section (@sections) {
		my ($pre, $desc);
		if (exists $INC{$path} && defined $INC{$path}) {
		    ($pre, $desc) = pod_section_to_array($section,$INC{$path});
		} elsif ( -f $plugin ) {
		    ($pre, $desc) = pod_section_to_array($section,$plugin);
		}
		$desc = $pre if (!defined $desc || $desc =~ /^\s*$/);
		push(@cols, $desc);
	    }

	    push(@desc, [ basename($name), @cols, $plugin ]);
	}
    }

    # Nothing more to do...
    return @desc;
}

=head2 pod_section_to_array

 Title   : pod_section_to_array
 Usage   : my ($name, $desc) = pod_section_to_array($section, $path)
 Function: Extract a section of a POD file as text
 Returns : string
 Args    : section name, POD file path

=cut

sub pod_section_to_array {
    my ($section,$path) = @_;

    my $name = undef;
    my $description = "";
    my $in_section = 0;
    open(my $pod, "<$path");
    while (<$pod>) {
	chomp;
	last if ($in_section && /^\=(head|cut|for)/);
	if (/^=head1\s*${section}\b/) {
	    $in_section = 1;
	    next ;
	}
	next if (!$in_section || /^\s*$/);
	s/^\s*//; s/\s*$//; # Remove trailing spaces
	$description .= length $description ? " $_" : $_;
    }
    close($pod);
    ($name, $description) = split(/\s*\-+\s*/,$description,2);

    return ($name, $description);
}

=head2 ruler

 Title   : ruler
 Usage   : ruler($start,$end,$marks,$interval,$reverse)
 Function: dump data structure as text
 Returns : string
 Args    : four integers and one boolean

           FIRST and SECOND arguments are the start
           and end coordinates of the interval you
           want the scale to represent

           The THIRD argument is the length of the 
           subdivisions or intervals. At the end of
           each interval, a mark (a colon, ":") is
           printed

           Marks will be replaced by numbers whenever
           coordinates are multiples of the FOURTH 
           argument

           The FIFTH argument indicates whether or not
           the ruler should be reversed (end to start)

=cut

sub ruler {
    my $self = shift if (ref($_[0]) eq __PACKAGE__ || $_[0] eq __PACKAGE__);
    my ($start, $end, $marks, $number, $reverse) = @_;

    # Initialize
    my $length = $end-$start+1;
    my $pad    = $marks-$start%$marks;
    my $ruler  = sprintf("%-${pad}s",$start);
    $ruler     = sprintf("%${pad}s",$start) if (defined $reverse && $reverse);

    # Add intervals
    $pad = $marks;
    foreach my $pos (int($start/$marks)+1..int($end/$marks)) {
	my $i = ($pos*$marks%$number == 0 ? $pos*$marks : ":");
	$pad = $length-length($ruler) if (length($ruler)+$pad > $length);
	$i = "|" if ($pad && length($i) > $pad);
	if (defined $reverse && $reverse) {
	    $ruler = sprintf("%${pad}s",$i) . $ruler;
	} else {
	    $ruler .= sprintf("%-${pad}s",$i);
	}
	$pad = $marks;
    }

    return $ruler;
}

=head2 ref2text

 Title   : ref2text
 Usage   : ref2text($ref)
 Function: dump data structure as text
 Returns : hash/array reference
 Args    : none

=cut

sub ref2text {
    my ($ref) = @_;
    return "<undef>" unless (defined $ref);
    if (ref $ref) {
	if (ref $ref eq "HASH") {
	    $ref = "{ ".join(", ",map { "$_ => ".ref2text($ref->{$_}) } keys %$ref)." }";
	} elsif (ref $ref eq "ARRAY") {
	    $ref = "[ ".join(", ",map { ref2text($_) } @$ref)." ]";
	} elsif ($ref eq "CODE") {
	} elsif ($ref eq "SCALAR") {
	    ref2text($ref);
	}
    }
    return $ref;
}

=head2 relative2absolute

 Title   : relative2absolute
 Usage   : $full_path = relative2absolute("data.txt")
 Function: get the absolute (fully qualified) path of a file/directory
 Returns : path to the file
 Args    : path to the file

=cut

sub relative2absolute {
    my ($path) = @_;
    $path = getcwd()."/$path" unless ($path =~ /^\//);
    while ($path =~ m|/\.{1,2}\/|) {
	$path =~ s|[^\/]+/\.\./||;
	$path =~ s|/\./||;
    }
    return $path;
}

=head2 nr2ids

 Title   : nr2ids
 Usage   : nr2ids($header)
 Function: parse FASTA header from NCBI's
           non-redundant databases 
 Returns : array of hashes
 Args    : full FASTA header (string)

Returned hash keys: see _new_dbxref in this code

=cut

sub nr2ids {
    my $self = shift if (ref($_[0]) eq __PACKAGE__ || $_[0] eq __PACKAGE__);
    my $opts = shift if (defined $_[0] && UNIVERSAL::isa($_[0],"HASH"));
    my ($header) = shift;
    chomp($header);

    # Fix variations in non-redundant input: make header delimiters look like 'fastacmd -c T' output
    $header =~ s/^\>+//; # Make sure FASTA header marker is removed
    if ($header !~ /\cA/) {
	$header =~ s/ +\; +gi\|/\cAgi\|/g;             # Replace faunique delimiter with Control-A
	$header =~ s/\s*\>gi\|/\cAgi\|/g;              # Replace fastacmd and blastdbcmd's default ' >' delimiter
	$header =~ s/^\cA+//;
    }

    # Process each match
    my @aoh  = ();
    my %seen = ();
    foreach my $entry (split(/\cA/,$header)) {
	my ($id,$desc) = split(/\s+/,$entry,2);
	my $organism = undef;
	next if (exists $seen{$id});
	$seen{$id} = 1;
	my ($start,$end) = ($id =~ /[:\/\|](\d+)[\-\.]{1,2}(\d+)$/);
	$id =~ s/[:\/\|](\d+)[\-\.]{1,2}(\d+)$//;
	my @id = map { s/^\s+//; s/\s+$//; $_ } split(/\|/,$id);
	next unless (scalar @id);
	if (defined $desc) {
#	    $desc =~ s/\s*\>.*$//;
	    ($organism) = ($desc =~ /\s*\[(.+)\]\s*$/);
	    $desc =~ s/\s*\[.*\]\s*$//;
	}
	my $type = undef;

	# Solo
	if (scalar @id == 1) {
	    my $dbname = undef;
	    if ($id[0] =~ /^\d+$/) {
		($dbname,$type) = ("Genbank","GI");
	    } elsif ($id[0] =~ /[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}\.\d*/) {
		($dbname,$type) = ("UniProtKB","ACCESSION");
	    } elsif ($id[0] =~ /^[A-Z]+\d+\.\d+$/) {
		($dbname,$type) = ("Genbank","ACCESSION");
	    } elsif ($id[0] =~ /^[A-Z]+_\d+\.\d+$/) {
		($dbname,$type) = ("RefSeq","ACCESSION");
	    }
	    my $hash = _new_dbxref($id[0], $dbname, 0, $id[0], $desc, $type, $organism, $start, $end);
	    if ($hash->{accession} =~ /[A-Z]{2}_\d+\.?\d*/) {
		$hash->{dbname} = "RefSeq";
	    } elsif ($hash->{accession} =~ /[A-Z]{3}\d+\.?\d*/) {
		$hash->{dbname} = "GenBank";
	    }
	    $hash->{version} = ($hash->{accession} =~ /\.(\d+)$/) || 0;
	    push (@aoh, $hash);
	    next;
	}

	# First
	next unless (defined $id[0] && length $id[0]);
	next unless (defined $id[1] && length $id[1]);
	$id[0] = uc($id[0]);
	if ($id[1] =~ /^\d+$/) { # Looks like a kosher GI
	    ($id[0],$type) = $id[0] eq 'GI' ? ("Genbank","GI") : ($id[0],"GI");
	    push (@aoh, _new_dbxref($id[1], $id[0], 0, $id[3] || $id[1], $desc, $type, $organism, $start, $end));
	} else { # Be aware of fake GIs!!!!
	    my $version = 0;
	    if (exists $Sources{$id[0]} && $id[0] ne "GI") {
		($id[0],$type) = ($Sources{$id[0]},"ACCESSION");
	    } elsif ($id[1] =~ /^[A-Z][a-z]+\d+$/) {
		($id[0],$type) = ("FADB","fakeGI");
	    } else { # I don't know this one...
		$type = undef;
	    }
	    push(@aoh, _new_dbxref($id[1], $id[0], $version, $id[3] || $id[1], $desc, $type, $organism, $start, $end));
	    last if (!defined $type  || $type eq 'GI_fragment');
	    last if (!defined $id[0] || $id[0] eq 'FADB'); # Additional fields in FADB entries are too variable to be meaningful...
	}

	# Second
	next unless (defined $id[2] && length $id[2]);
	$id[2] = exists $Sources{$id[2]} ? $Sources{$id[2]} : next; # Stop if accession database is unknown
	$aoh[$#aoh]->{"dbname"} = $id[2];

	# Third
	$type = "ACCESSION";
	if (defined $id[3] && length $id[3]) {
	    my ($version) = ($id[3] =~ /\.(\d+)$/);
	    $version = 0 unless (defined $version && length $version);
	    if ($id[2] eq 'PDB') { # PDB
		$id[3] .= defined $id[4] ? "_$id[4]" : "_A"; # add chain
		push(@aoh, _new_dbxref($id[3], $id[2], $version, $id[3], $desc, $type, $organism, $start, $end));
	    } elsif ($id[2] eq 'SwissProt') {
		push(@aoh, _new_dbxref($id[3], $id[2], $version, $id[3], $desc, $type, $organism, $start, $end));
	    } else {
		push(@aoh, _new_dbxref($id[3], $id[2], $version, $id[3], $desc, $type, $organism, $start, $end));
	    }
	}

	# Fourth
	$type = 'LOCUS';
	if (defined $id[4] && length $id[4]) {
	    #next if ($id[2] eq 'PIR'); # PIR locus name may collide with NCBI nucleotide entries
	    next if ($id[2] eq 'PDB');
	    my ($version) = ($id[4] =~ /\.(\d+)$/);
	    push(@aoh, _new_dbxref($id[4], $id[2], defined $version ? $version : 0, $id[1], $desc, $type, $organism, $start, $end));
	}
    } # foreach my $id (@name)

    return @aoh;
}

# Subroutine to make sure all hashes look the same
sub _new_dbxref {
    return { 
	accession   => $_[0],
	dbname      => $_[1] || "UnknownSeqDB",
	version     => $_[2] || 0,
	accgroup    => $_[3] || "NO_GROUP",
	description => $_[4] || "",
	acctype     => $_[5] || "UnknownID",
	organism    => $_[6] || "NO_NAME",
	start       => $_[7] || 0,
	end         => $_[8] || 0,
    };
}

=head2 ids2nr

 Title   : ids2nr
 Usage   : ids2nr($ref)
 Function: convert parsed FASTA identifier
           to sequence identifiers
 Returns : array of strings
 Args    : array of hashes (see nr2ids)

 The first hash may be used to set the concatenation
 string for the FASTA headers (key 'nr', default: " >")

=cut

sub ids2nr {
    my $self = shift if (ref($_[0]) eq __PACKAGE__ || $_[0] eq __PACKAGE__);
    my $opts = shift if (defined $_[0] && !exists $_[0]->{accgroup});
    my $sep  = exists $opts->{concatenate} ? $opts->{concatenate} : " >";
    $opts->{description} = 1 unless (exists $opts->{description});

    # Group by accgroup
    my %group = ();
    foreach my $hash (@_) {
	push(@{ $group{$hash->{accgroup}} }, $hash);
    }

    # Order by: group, acctype, accession, version
    my @header = ();
    foreach my $groupid (keys %group) {
	my @hash = sort {
	    my $typeA = exists $SeqIdTypeRank{$a->{acctype}} ? $SeqIdTypeRank{$a->{acctype}} : scalar(@SeqIdTypes);
	    my $typeB = exists $SeqIdTypeRank{$b->{acctype}} ? $SeqIdTypeRank{$b->{acctype}} : scalar(@SeqIdTypes);
	    $typeA <=> $typeB || $b->{accession} <=> $a->{accession} || 
		$b->{accession} cmp $a->{accession} || $a->{version} <=> $b->{version}
	} @{ $group{$groupid} };

	# Collapse by group
	my @ids    = ();
	foreach my $hash (@hash) {
	    # GI
	    if (grep { $_ eq $hash->{acctype} } qw(GI fakeGI GI_fragment)) {
		push(@ids, 'gi');
	    }

	    # Accession
	    if ($hash->{acctype} eq 'ACCESSION') {
		if (exists $ReversedSources{$hash->{dbname}}) {
		    push(@ids, lc $ReversedSources{$hash->{dbname}});
		} else {
		    carp "Unknown dbname ".$hash->{dbname};
		    push(@ids, lc($hash->{dbname}));
		}
	    }

	    if ($hash->{dbname} eq 'PDB' && !looks_like_number($hash->{accession})) {
		my ($chain) = $hash->{accession} =~ /([A-Z]+)$/;
		$hash->{accession} =~ s/${chain}$/\|${chain}/;
		$hash->{accession} = uc $hash->{accession};
	    }

	    push(@ids, $hash->{acctype} eq 'ACCESSION' && $hash->{version} ? $hash->{accession}.".".$hash->{version} : $hash->{accession});
	}

	my $text = join("|",@ids);
	$text .= "|" unless ($hash[$#hash]->{acctype} eq "LOCUS" || $hash[0]->{dbname} eq "PDB");
	$text .= " ".$hash[0]->{description} if ($opts->{description} && exists $hash[0]->{description});
	push(@header, $text);
    }

    my $header = join($sep, @header);
    return $header;
}

# Make perl compiler happy
1;
