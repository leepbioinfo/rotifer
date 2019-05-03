=head1 NAME

Rotifer::Parser::genbank2table - load Genbank data into a sparse table using Bio::SeqIO

=head2 SYNOPSIS

 my $parser = Rotifer::Parser->create(format => "genbank2table");
 foreach my $table ($parser->next) {
   ... do something ...
 }

=head2 DESCRIPTION

Rotifer::Parser::genbank2table is a parser based on Bioperl
(Bio::SeqIO::RichSeqI) that loads sequence annotations into
a table (Rotifer::Data::Table).

The parse uses a iterator routine to data extraction but the
user has to request the type of information that should be
extracted. The intention is to 

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Rotifer::Parser::genbank2table;

use Data::Dumper;
use Bio::SeqIO;
use Moose;
use strict;
use warnings;

=head2 ATTRIBUTES

=head2 annotation

 Name    : annotation
 Function: tags for whole sequence comments, references or
           other annotations
 Returns : array reference
 Args    : array reference

=cut

has 'annotation' => (
    is      => 'rw',
    isa     => 'ArrayRef[Str]',
    default => sub { [] },
);

=head2 clean

 Name    : clean
 Function: List of regular expression to clean column names. 
           The regular expressions will also be parsed and used
           to name the correponding columns if matching /(\S+):$/
 Returns : array reference
 Args    : array reference

=cut

has 'clean' => (
    is      => 'rw',
    isa     => 'ArrayRef[Str]',
    default => sub { [] },
);

=head2 columns

 Name    : columns
 Function: annotation tags to extract
 Returns : array reference
 Args    : array reference

 Note    : you may append regular expression to narrow
           the tag selection procedure. E.g.

             columns => [ 'db_xref' ]

           will load all db_xref entries but

             columns => [ 'db_xref=GI:' ]

           will extract only db_xref values matching the GI:
           regular expression.

           Using 

             columns => [ 'all' ]

           will load all tags. On the other hand,

             columns => [ 'all', 'db_xref=GI:' ]

           will load all tags but keep GI values in a separate
           column (i.e. hash key) but note that the order the
           columns are specified matters: more specific regular
           expressions should be given first!!!!

=cut

has 'columns' => (
    is      => 'rw',
    isa     => 'ArrayRef[Str]',
    default => sub { [ 'product', 'db_xref=GI:', 'protein_id', 'db_xref=GeneID:', 'gene', 'locus_tag' ] },
);

=head2 exclude

 Name    : exclude
 Function: list of features that should be ignored
 Returns : array reference
 Args    : array reference

=cut

has 'exclude' => (
    is      => 'rw',
    isa     => 'ArrayRef[Str]',
    default => sub { [] },
);

=head2 format

 Name    : format
 Function: Input file format (see Bio::SeqIO)
 Returns : array reference
 Args    : array reference

=cut

has 'format' => (
    is      => 'rw',
    isa     => 'Str',
    default => 'genbank',
);

=head2 include

 Name    : include
 Function: list of features to extract
 Returns : array reference
 Args    : array reference

=cut

has 'include' => (
    is      => 'rw',
    isa     => 'ArrayRef[Str]',
    default => sub { [] },
);

=head2 loaded_files

 Name    : loaded_files
 Function: list of input files loaded
 Returns : array reference
 Args    : array reference

 Handlers:

  all_loaded_files : retrieve all file names
  add_loaded_files : add one or more values
  nof_loaded_files : number of files loaded

=cut

has 'loaded_files' => (
    traits  => ['Array'],
    is      => 'rw',
    isa     => 'ArrayRef[Str]',
    default => sub { [] },
    handles => {
       all_loaded_files => 'elements',
       add_loaded_files => 'push',
       nof_loaded_files => 'count',
    }
);

=head2 preferred_taxons

 Name    : preferred_taxons
 Function: List of taxon names to use in a compact representation
           of organismal taxonomy
 Returns : array reference
 Args    : array reference

=cut

has 'preferred_taxons' => (
    is      => 'rw',
    isa     => 'ArrayRef[Str]',
    default => sub { [] },
);

=head2 remove_columns

 Name    : remove_columns
 Function: list of tags that will NOT be extracted
 Returns : array reference
 Args    : array reference

=cut

has 'remove_columns' => (
    is      => 'rw',
    isa     => 'ArrayRef[Str]',
    default => sub { [] },
);

=head2 seqio_options

 Usage   : Rotifer::Parser->new(seqio_options => { "-file" => "a.fa" })
 Function: get/set Bio::SeqIO options
 Returns : hash reference
 Value   : hash
 Builder : _default_seqio_options
 Trigger : _seqio_options_trigger

=cut

has 'seqio_options' => (
    is      => "rw", 
    isa     => "HashRef",
    lazy    => 1,
    builder => '_default_seqio_options',
    trigger => \&_seqio_options_trigger,
    );

sub _default_seqio_options {
    return {};
}

sub _seqio_options_trigger {
    my ($self, $new, $old) = @_;
    foreach my $key (keys %$new) {
	next if ($key =~ /^-+/);
	my $value = delete $new->{"$key"};
	$new->{"-$key"} = $value;
    }
}

=head2 METHODS

=head2 parse

 Title   : parse
 Usage   : $parser->parse(\%opts, @files)
 Function: parse genbank to hashes
 Returns : two array references and a list of hash references

           The array references are the name of the columns
           output columns and the keys used to find 
           corresponding values in the hashes

 Args    : hash

=cut

sub parse {
    my $self = shift;

    # Pre-process options
    my %preferred      = map { ($_,1) } @{$self->preferred_taxons};
    my %clean          = map { ($_,1) } @{$self->clean};
    my %exclude        = map { ($_,1) } @{$self->exclude};
    my %remove_columns = map { ($_,1) } @{$self->remove_columns};
    my %annotation     = map { ($_,1) } @{$self->annotation};

    # Parse and index tags
    my %tags = ();
    my @columns = ();
    foreach my $tag (@{ $self->columns }) {
	next if (exists $remove_columns{$tag});
	my ($name,$patt) = split(/=/,$tag,2);
	next if (exists $remove_columns{$name});
	my $has_regexp = defined $patt || $name eq 'all';
	for (1..2) {
	    if ($_ == 2) {
		last if ($has_regexp);
		$patt = "${name}:";
		$name = 'db_xref';
	    } else {
		$patt = ".+" if (!defined $patt);
	    }
	    if (exists $tags{$name}) {
		push(@{ $tags{$name} }, $patt);
	    } else {
		$tags{$name} = [ $patt ];
	    }
	    push(@columns, $name);
	}
    }
    unshift(@{$self->include},"source") 
	if (scalar @{$self->annotation} && !grep { $_ eq "source" } @columns);

    # Load and parse files
    my %keys  = ();
    my @table = ();
    my @base  = qw(accession primary_id);
    my @fbase = qw(primary_tag start end strand);
    foreach my $file (@_) {
	# Open input file
	my @file = ();
	if (ref $file eq "GLOB") {
	    @file = (-fh => $file);
	} elsif (! -f $file && $file eq "-") {
	    $self->add_loaded_files($file);
	    @file = (-fh => \*STDIN);
	} else {
	    $self->add_loaded_files($file);
	    @file = (-file => "<$file");
	}
	my $in = Bio::SeqIO->new(-format=>$self->format, @file, %{$self->seqio_options});

	# Process file
	# Set tags
	while (my $seq = $in->next_seq) {
	    my $hash = { map { ($_, $seq->$_()) } @base };
	    $hash->{locus} = $seq->display_id;

	    # Sequence you want annotations...
	    if (scalar @{$self->annotation}) {
		# Sequence properties
		map { $hash->{$_} = $seq->$_() } qw(molecule alphabet description);

		# Taxonomy
		if (defined $seq->species) {
		    my %gettag = (
			ncbi_taxid      => 'ncbi_taxon_id',
			scientific_name => 'name',
			);
		    map {
			my $tag = $_;
			$tag = exists $gettag{$_} ? $gettag{$_} : $_;
			$hash->{$tag}  = $seq->species->taxon->$_();
		    } keys %gettag;
		    my @lineage = $seq->species->classification;
		    if (scalar @{$self->preferred_taxons}) {
			@lineage = grep { exists $preferred{$_} } map { lc($_) } @lineage;
			$hash->{lineage}  = join(">",@lineage);
		    } else {
			$hash->{classification} = join(", ",@lineage);
		    }
		}

		# Annotation
		my %count = ();
		foreach my $ann ($seq->annotation->flatten_Annotations) {
		    my $tag = $ann->tagname;
		    next unless (exists $annotation{all} || exists $annotation{$tag});

		    # Retrieve value
		    my $value = undef;
		    if ($ann->isa("Bio::Annotation::SimpleValue")) {
			$value = $ann->value;
		    } elsif ($ann->isa("Bio::Annotation::StructuredValue")) {
			$value = join(", ",$ann->get_all_values);
		    } elsif ($ann->isa("Bio::Annotation::Comment")) {
			$value = $ann->value;
		    } elsif ($ann->isa("Bio::Annotation::Reference")) {
			$tag   = "db_xref\cAPMID";
			$value = $ann->pubmed if (defined $ann->pubmed);
		    } elsif ($ann->isa("Bio::Annotation::DBLink")) {
			$value = defined $ann->version ? $ann->primary_id.".".$ann->version : $ann->primary_id;
			next if (!defined $ann->database);
			my $db = $ann->database;
			if ($db eq "PDB") {
			    my @val = split(/\_/,$value);
			    $val[0] =~ tr/A-Z/a-z/;
			    $value = join("",@val);
			} elsif ($ann->database eq "GenPept") {
			    $db = $value =~ /^[A-Z][A-Z]_/ ? 'RefSeq' : 'Genbank';
			}
			$tag = "db_xref\cA".$db;
		    }

		    # If there are values
		    if (defined $value && length $value) {
			$value =~ s/\.$//;
			if (exists $hash->{$tag}) {
			    $tag .= "\cA".$count{$tag}++;
			    $hash->{$tag} = $value;
			} else {
			    $hash->{$tag} = $value;
			}
		    }
		} # foreach my $ann ($seq->annotation->flatten_Annotations)
	    } # if (scalar @{$self->annotation})

	    # Features
	    my $has_features = 0;
	  FEAT: foreach my $generic ($seq->get_all_SeqFeatures) {
	      # Check whether this feature was requested
	      next if (exists $exclude{$generic->primary_tag});
	      unless (grep { $_ eq 'all' } @{$self->include}) {
		  my $matches = 0;
		  foreach my $tag (@{$self->include}) {
		      $matches++ if ($generic->primary_tag =~ /$tag/);
		  }
		  next unless ($matches);
	      }
#	      next if (scalar @{$self->annotation} && $generic->primary_tag ne "source");

	      # Initialize data for this feature: tag - value annotations
	      my $feat = { %$hash, map { ($_, $generic->$_()) } @fbase };

	      if (scalar @columns) {
		  # Tag-based interface
		  my @tags = exists $tags{all} ? $generic->get_all_tags : @columns;
		  foreach my $tag (@tags) {
		      next if (exists $remove_columns{$tag});
		      next unless $generic->has_tag($tag);
		      my @regexps = @{$tags{$tag}} if (exists $tags{$tag});
		      push(@regexps,'.+') if (exists $tags{all});
		      my %loaded = ();
		      foreach my $regexp (@regexps) {
			  next FEAT if (exists $exclude{$generic->primary_tag . "=" . $tag . "=" . $regexp});
			  next FEAT if (exists $exclude{$generic->primary_tag . "=" . $tag . "=" . $regexp . ":"});
			  my @values = sort grep { /$regexp/ } $generic->each_tag_value($tag);
			  @values = grep { !exists $loaded{"$tag\cA$_"} } @values;
			  next unless (scalar @values);
			  map { $loaded{"$tag\cA$_"} = 1 } @values;
			  my $key = $tag;
			  if ($regexp ne '.*' && $regexp ne '.+') {
			      @values = map { $loaded{"$key\cA$_"} = 1; s/$regexp//; $_ } @values if (exists $clean{$tag});
			      $key = "$tag\cA$regexp";
			  }
			  $feat->{$key} = join(" | ",map { s/_no_value/1/; $_ } @values) if (grep { length } @values); # Feature has tag matching $tags{$tag}
		      }
		  }

		  # Annotation collection-based interface
		  if ($generic->can("annotation")) {
		      my @anntags = exists $tags{all} ? $generic->annotation->get_all_annotation_keys : @columns;
		      foreach my $tag (@anntags) {
			  next if (exists $remove_columns{$tag});
			  my @ann = $generic->annotation->flatten_Annotations($tag);
			  next unless (scalar @ann);
			  my @regexps = @{$tags{$tag}} if (exists $tags{$tag});
			  push(@regexps,'.+') if (exists $tags{all});
			  foreach my $regexp (@regexps) {
			      next FEAT if (exists $exclude{$generic->primary_tag . "=" . $tag . "=" . $regexp});
			      next FEAT if (exists $exclude{$generic->primary_tag . "=" . $tag . "=" . $regexp . ":"});
			      my %loaded = ();
			      foreach my $ann (@ann) {
				  my @values = sort grep { /$regexp/ } map { $_->{'value'} } $ann->hash_tree;
				  @values = grep { !exists $loaded{"$tag\cA$_"} } @values;
				  next unless (scalar @values);
				  map { $loaded{"$tag\cA$_"} = 1 } @values;
				  my $key = $tag;
				  if ($regexp ne '.*' && $regexp ne '.+') {
				      @values = map { s/$regexp//; $_ } @values if (exists $clean{$tag});
				      $key = "$tag\cA$regexp";
				  }
				  $feat->{$key} = join(" | ",map { s/_no_value/1/; $_ } @values); # Feature has tag matching $tags{$tag}
			      }
			  }
		      }
		  }
	      } # if (scalar @columns)

	      # Register output columns
	      map { $keys{$_} = 1 } keys %{$feat};
	      push(@table, $feat);
	      $has_features = 1;
	  } # FEAT: foreach my $generic ($seq->get_all_SeqFeatures)

	    push(@table, $hash) unless ($has_features);
	} # while (my $seq = $in->next_seq)

	# Close input
	$in->close; undef($in);
    } # foreach my $file (@ARGV)

    # Prepare header
    push(@base, 'locus');
    my %map = (
	primary_id  => 'seqid',
	primary_tag => 'feature',
	);
    my @header = map { delete $keys{$_}; $map{$_} || $_ } (@base,@fbase);
    push(@header, 
	 map { 
	     my ($tag,$re) = split(/\cA/);
	     if (defined $re && $re ne '.*' && $re ne '.+' && $re ne '.+:') {
		 $re =~ s/:$//;
		 exists $map{$re} ? $map{$re} : $re =~ /^\d+$/ ? $tag : $re;
	     } else {
		 exists $map{$tag} ? $map{$tag} : $tag;
	     }
	 } sort keys %keys);
    my @keys = (@base, @fbase, sort keys %keys);

    return (\@header, \@keys, @table);
}

__PACKAGE__->meta->make_immutable;
1;
