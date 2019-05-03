# POD documentation - main docs before the code

=head1 NAME

Rotifer::DB::FADB - method to handle sequence databases
                    using site-specific identifiers.

=head1 SYNOPSIS

  # Retrieving FASTA sequences
  use Rotifer::DB::FADB qw(fakegi2fasta);
  my $fasta = fakegi2fasta($fakegi);

=head1 DESCRIPTION

Rotifer::DB::FADB provides a simple interfaces to retrieve
data from local FASTA, BLAST and other types of databases. 

=head1 CONFIGURATION

Configuration parameters are loaded from

 ~/.rotifer/etc/db/fadb.yml 

or, the file above doesn't exist, 

 <rotifer installation root>/etc/rotifer/db/fadb.yml

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::DB::FADB;

use autodie qw(:all);
use Carp::Clan qr/^Rotifer::DB::FADB/;
use File::Monitored qw(:temp);
use List::MoreUtils qw(uniq);
use Rotifer::DB::Base;
use Scalar::Util qw(reftype);
use Sub::Exporter -setup => {
    exports => [qw/fakegi2fasta is_fakeGI/],
};
use strict;
use warnings;
use vars qw(@FAKEGI_ORGANISMS $FAKEGI_ORGS_RE %FAKEGI_FILES $fakeGIs);

# Global variables
BEGIN {
#    $fakeGIs = "$ENV{ROTIFER_DATA}/taxonomy/fakegi_names.txt";
    $fakeGIs = Rotifer::DB::Base->data_path("taxonomy","fakegi_names.txt");
};

=head2 fakegi2fasta

 Title   : fakegi2fasta
 Usage   : Rotifer::DB->fakegi2fasta(\%opts, @fakegis)
 Function: use grep to retrieve sequences foram FASTA files
 Returns : file name or file handle
 Args    : (optional) hash reference and 
           array of sequence identifiers (strings)

=cut

sub fakegi2fasta {
    my $self = shift if (defined $_[0] && ($_[0] eq __PACKAGE__ || UNIVERSAL::isa($_[0],__PACKAGE__)));
    my $opts = { %{shift(@_)} } if (defined $_[0] && UNIVERSAL::isa($_[0],"HASH"));
    my @ids  = @_; @ids = map { s/^\s+//; s/\s+$//; $_ } @ids;

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

    # Classify IDs by organism abbreviation
    my %abbr = ();
    foreach my $id (@ids) {
	my ($abbr) = ($id =~ /^([A-Za-z]+)\d+/);
	push(@{$abbr{$abbr}}, $id) if (defined $abbr);
    }

    # Grep each fakeGI organism file
    my @fasta = ();
    foreach my $abbr (keys %abbr) {
	my ($fh,$temp) = &ttempfile();
	print $fh join("\n",@{$abbr{$abbr}}),"\n";
	close($fh);
#	my $file = "$ENV{ROTIFER_DATA}/fadb/euk_fa/fakegi/$FAKEGI_FILES{$abbr}";
	my $file = Rotifer::DB::Base->data_path("fadb","euk_fa","fakegi",$FAKEGI_FILES{$abbr});
	die "No such fake GI FASTA file $file" unless (defined $file);
	open(my $grep, "grep -A 1 -w -F -f $temp $file |" );
	push(@fasta,grep { !/^[\-\s]*$/ } <$grep>);
	print $outfh join('',@fasta);
	close($grep);
	unlink($temp);
    }

    return $output;
}

=head2 is_fakeGI

 Title   : is_fakeGI
 Usage   : Rotifer::DB::FADB->is_fakeGI($id)
 Function: see if input looks like a fake GI
 Returns : boolean
 Args    : string

=cut

sub is_fakeGI {
    my $self = shift if (defined $_[0] && ($_[0] eq __PACKAGE__ || UNIVERSAL::isa($_[0],__PACKAGE__)));
    my $id = shift;
    return 0 unless (defined $id);

    if (!defined $FAKEGI_ORGS_RE) {
	my @fakeorgs = &fakeGI_organisms;
	return 0 unless (scalar @fakeorgs);
	$FAKEGI_ORGS_RE = join("|",sort { lc($a) cmp lc($b) } map { $_->{abbreviation} } @fakeorgs);
	$FAKEGI_ORGS_RE = qr/^($FAKEGI_ORGS_RE)/o;
    }

    return ($id =~ m/$FAKEGI_ORGS_RE/) ? 1 : 0;
}

=head2 fakeGI_organisms

 Title   : fakeGI_organisms
 Usage   : Rotifer::DB::FADB->fakeGI_organisms(@files)
 Function: Retrieve organisms whose sequences
           are identified by fake GIs
 Returns : array of hash references
 Args    : one or more table file names

=cut

sub fakeGI_organisms {
    my $self = shift if (defined $_[0] && ($_[0] eq __PACKAGE__ || UNIVERSAL::isa($_[0],__PACKAGE__)));

    if (!scalar @FAKEGI_ORGANISMS) {
	if (! -f $fakeGIs) {
	    carp "Unable to read fakeGIs file $fakeGIs";
	    return ();
	}

	# Parse fakegi organisms data
	my @columns = qw(abbreviation file name preferred lineage);
	foreach my $file ($fakeGIs) {
	    open(my $fh,"<$file");
	    while (<$fh>) {
		chomp;
		my @data = map { s/^\s+//; s/\s+$//; $_ } split(/\s+:\s+/);
		my $hash =  { map { ($columns[$_],$data[$_]) } (0..3,3) };
		push(@FAKEGI_ORGANISMS,$hash);
		$FAKEGI_FILES{$hash->{abbreviation}} = $hash->{file};
	    }
	    close($fh);
	}
    }

    return @FAKEGI_ORGANISMS;
}

# Make perl compiler happy
1;
