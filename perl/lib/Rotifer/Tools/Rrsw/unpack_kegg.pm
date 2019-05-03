# POD documentation - main docs before the code

=head1 NAME

Rotifer::Tools::Rrsw::unpack_kegg - unpack, split and format Kegg databases

=head1 SYNOPSIS

  # Creating a new tool script
  use Rotifer::Tools::Rrsw;
  my $rrsw = Rotifer::Tools::Rrsw->create("unpack_kegg", { name => $target });
  $rrsw->get(@ARGV);

=head1 DESCRIPTION

Rotifer::Tools::Rrsw::unpack_kegg unpacks Kegg database files and prepares
these files to be used in searches with HMMER 3 and BLAST.

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::Tools::Rrsw::unpack_kegg;

use strict;
use warnings;
use autodie qw(:all);
use Carp::Clan qr/^Rotifer::Tools::Rrsw/;
use File::Basename;
use File::Which;
use IPC::Run qw(run);
use Moose;
use Rotifer;
extends 'Rotifer::Tools::Rrsw::unpack';

=head2 target

 Usage   : Rotifer::Tools::Rrsw::unpack_kegg->new(target => $target)
 Function: basename of the target
 Value   : string
 Required: 1

=cut

# Inherited from Rotifer::Tools::Rrsw::unpack

=head2 test

 Usage   : Rotifer::Tools::Rrsw->new(test => 1)
 Function: run in test mode
 Value   : Boolean

=cut

# Inherited from Rotifer::Tools::Rrsw::unpack

=head2 options

 Usage   : Rotifer::Tools::Rrsw::unpack_kegg->new(options => { gzip => "/usr/bin/gzip" })
 Function: module configuration
 Value   : hash reference as returned by Config::General
 Builder : _default_options

=cut

# Inherited from Rotifer::Tools::Rrsw::unpack

=head2 executable

 Usage   : Rotifer::Tools::Rrsw::unpack_kegg->new(executable => { gzip => "/usr/bin/gzip" })
 Function: path to tar, gzip and other executable
 Value   : hash reference with program names as keys
 Builder : _default_executable (uses File::Which)

=cut

sub _default_executable {
    my $self = shift;
    my $paths = $self->next::method();
    my $path = which('makeblastdb');
    $paths->{makeblastdb} = $path if (defined $path);
    return $paths;
}

=head2 suffixes

 Usage   : Rotifer::Tools::Rrsw::unpack_kegg->new(suffixes => { '.tar.gz' => "tar" })
 Function: map of suffixes to unpacking programs
 Value   : hash reference with filename suffixes as keys
 Builder : _default_executable (uses File::Which)

=cut

sub _default_suffixes {
    my ($self) = @_;
    my $ext = $self->next::method;
    $ext->{'.pep'} = 'makeblastdb';
    return $ext;
}

=head2 METHODS

=head2 execute

 Title   : execute
 Usage   : $count = $rrsw->execute($source, $dest)
 Function: unpack downloaded files
 Returns : list of files unpacked
 Args    : input file(s)/directory and output directory

=cut

sub execute {
    my ($self, $source, $unpack) = @_;
    my $status = 0;

    # Unpack first
    $self->next::method($source, $unpack);

    # Detect appropriate executable
    my ($cmd,$ext) = $self->get_executable_by_extension("anything.pep");
    $cmd = $self->executable->{$cmd} if (defined $cmd && exists $self->executable->{$cmd});
    my @basecmd = ( $cmd, "-dbtype", "prot" );

    # Laod data about this release
    die "Could not cd to  unpack directory $unpack" if ( ! -d $unpack );
    chdir($unpack) || die "Could not cd to";
    open(RELEASE,"<RELEASE") || die "Could not read release file";
    chomp(my $release = <RELEASE>);
    close(RELEASE);

    # Laod organism identifiers
    my %tax = ();
    open(TAX,"<genes/taxonomy") || die "Could not read taxonomy file";
    while (<TAX>) {
	chomp;
	next if /^\#/;
	my @F = split(/\t/);
	$tax{$F[1]} = [ @F ];
    }
    close(TAX);

    # Process master FASTA file
    mkdir("genes/KEGG_fasta_database") if ( ! -d "genes/KEGG_fasta_database");
    chdir("genes/KEGG_fasta_database") || die "Failed to chdir to genes/KEGG_fasta_database directory";
    die "Could not find/read genes.pep" if ( ! -f "../fasta/genes.pep");
    unlink("all_orgs") if ( -l "all_orgs" );
    symlink("../fasta/genes.pep","all_orgs") || die "Could not symlink original genes.pep file";
    my @cmd = ( @basecmd, "-in", "all_orgs", "-out", "all_orgs", "-title", "$release", "-logfile", "all_orgs.log" );
    $self->_run(@cmd);

    # Split and format genes.pep
    my $organism = "";
    open(GENES,"<all_orgs");
    while (<GENES>) {
	if (my ($new) = /^>(\S+):/) {  # FASTA header found!
	    if ($organism ne $new) {   # Is it a new organism?
		if ($organism ne "") { # Is it the first organism?
		    close(ORG) if (!$self->test);
		    my $name = exists $tax{$organism} ? $tax{$organism}->[3] || $organism : $organism;
		    my @cmd = ( @basecmd, "-in", "$organism", "-out", "$organism", "-title", "${name} proteome in $release", "-logfile", "${organism}.log" );
		    $self->_run(@cmd);
		}
		if (!$self->test) {
		    open(ORG,">$new") || die "Could not create genome proteome file $new";
		}
	    } # if ($organism ne $new)
	    $organism = $new;
	}
	print ORG $_ if (!$self->test);
    }
    close(GENES);

    # Last organism!
    if ($organism ne "") {
	close(ORG) if (!$self->test);
	my $name = exists $tax{$organism} ? $tax{$organism}->[3] || $organism : $organism;
	my @cmd = ( @basecmd, "-in", "$organism", "-out", "$organism", "-title", "${name} proteome in $release", "-logfile", "${organism}.log" );
	$self->_run(@cmd);
    }

    # Go to initial directory and return status
    chdir($Rotifer::IWD);
    return $status;
}

sub _run {
    my ($self, @cmd) = @_;

    my $status = 0;
    if ($self->test) {
	print join(" ",@cmd),"\n";
    } else {
	$status = run(\@cmd);
    }

    return $status;
}

__PACKAGE__->meta->make_immutable;
1;
