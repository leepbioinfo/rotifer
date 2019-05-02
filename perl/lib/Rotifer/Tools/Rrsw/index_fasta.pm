# POD documentation - main docs before the code

=head1 NAME

Rotifer::Tools::Rrsw::index_fasta - unpack and format Blast databases

=head1 SYNOPSIS

  # Creating a new tool script
  use Rotifer::Tools::Rrsw;
  my $rrsw = Rotifer::Tools::Rrsw->create("index_fasta", { name => $target });
  $rrsw->get(@ARGV);

=head1 DESCRIPTION

Rotifer::Tools::Rrsw::index_fasta prepares FASTA files
to be used in searches by HMMER3 (using perl and esl-sfetch).

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::Tools::Rrsw::index_fasta;

use strict;
use warnings;
use autodie qw(:all);
use Carp::Clan qr/^Rotifer::Tools::Rrsw/;
use File::Basename;
use File::Find;
use File::Which;
use IPC::Run qw(run);
use Moose;
use Rotifer;
extends 'Rotifer::Tools::Rrsw::unpack';

=head2 target

 Usage   : Rotifer::Tools::Rrsw::index_fasta->new(target => $target)
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

 Usage   : Rotifer::Tools::Rrsw::index_fasta->new(options => { gzip => "/usr/bin/gzip" })
 Function: module configuration
 Value   : hash reference as returned by Config::General
 Builder : _default_options

=cut

# Inherited from Rotifer::Tools::Rrsw::unpack

=head2 executable

 Usage   : Rotifer::Tools::Rrsw::index_fasta->new(executable => { gzip => "/usr/bin/gzip" })
 Function: path to tar, gzip and other executable
 Value   : hash reference with program names as keys
 Builder : _default_executable (uses File::Which)

=cut

sub _default_executable {
    my $self = shift;
    my $paths = $self->next::method();
    my $path = which('blastdbcmd');
    $paths->{blastdbcmd} = $path if (defined $path);
    $path = which('esl-sfetch');
    $paths->{"esl-sfetch"} = $path if (defined $path);
    return $paths;
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

    # Unpack first
    $self->next::method($source, $unpack);

    # Format unpacked databases
    mkdir($unpack) if ( ! -d $unpack );
    chdir($unpack);

    # Detect appropriate executables
    my $eslsfetch = $self->executable->{'esl-sfetch'};
    die "Could not locate esl-sfetch. You should install HMMER with the configure flag --with-easel."
	unless (defined $eslsfetch);

    # Run
    my $status = 0;
    foreach my $db (@{$self->options->{'fasta'}}) {
	next if ($db =~ /\.\d\d$/); # Avoid fragments
	my @cmd = ([ $eslsfetch, "--index", "${db}" ]);
	if ($self->test) {
	    print join(" ",@{$cmd[0]}),"\n";
	} else {
	    $status = run(@cmd) if ( -f "${db}" );
	}
    }

    chdir($Rotifer::IWD);
    return $status;
}

__PACKAGE__->meta->make_immutable;
1;
