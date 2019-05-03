# POD documentation - main docs before the code

=head1 NAME

Rotifer::Tools::Rrsw::unpack_pfam - unpack and format PFam databases

=head1 SYNOPSIS

  # Creating a new tool script
  use Rotifer::Tools::Rrsw;
  my $rrsw = Rotifer::Tools::Rrsw->create("unpack_pfam", { name => $target });
  $rrsw->get(@ARGV);

=head1 DESCRIPTION

Rotifer::Tools::Rrsw::unpack_pfam unpacks PFam database files and prepares
these files to be used in searches by HMMER 3 (using hmmpress).

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::Tools::Rrsw::unpack_pfam;

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

 Usage   : Rotifer::Tools::Rrsw::unpack_pfam->new(target => $target)
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

 Usage   : Rotifer::Tools::Rrsw::unpack_pfam->new(options => { gzip => "/usr/bin/gzip" })
 Function: module configuration
 Value   : hash reference as returned by Config::General
 Builder : _default_options

=cut

# Inherited from Rotifer::Tools::Rrsw::unpack

=head2 executable

 Usage   : Rotifer::Tools::Rrsw::unpack_pfam->new(executable => { gzip => "/usr/bin/gzip" })
 Function: path to tar, gzip and other executable
 Value   : hash reference with program names as keys
 Builder : _default_executable (uses File::Which)

=cut

sub _default_executable {
    my $self = shift;
    my $paths = $self->next::method();
    my $path = which('hmmpress');
    $paths->{hmmpress} = $path if (defined $path);
    return $paths;
}

=head2 suffixes

 Usage   : Rotifer::Tools::Rrsw::unpack_pfam->new(suffixes => { '.tar.gz' => "tar" })
 Function: map of suffixes to unpacking programs
 Value   : hash reference with filename suffixes as keys
 Builder : _default_executable (uses File::Which)

=cut

sub _default_suffixes {
    my ($self) = @_;
    my $ext = $self->next::method;
    $ext->{'.hmm'} = 'hmmpress';
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

    # Unpack first
    $self->next::method($source, $unpack);

    # Format unpacked databases
    mkdir($unpack) if ( ! -d $unpack );
    chdir($unpack);

    # Detect appropriate executable
    my ($cmd,$ext) = $self->get_executable_by_extension("anything.hmm");
    $cmd = $self->executable->{$cmd} if (defined $cmd && exists $self->executable->{$cmd});

    # Run
    my $status = 0;
    foreach my $hmm (glob('*.hmm')) {
	my @cmd = ([ $cmd, "-f", $hmm ]);
	if ($self->test) {
	    print join(" ",@{$cmd[0]}),"\n";
	} else {
	    $status = run(@cmd);
	}
    }

    chdir($Rotifer::IWD);
    return $status;
}

__PACKAGE__->meta->make_immutable;
1;
