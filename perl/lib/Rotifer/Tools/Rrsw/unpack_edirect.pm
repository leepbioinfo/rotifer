# POD documentation - main docs before the code

=head1 NAME

Rotifer::Tools::Rrsw::unpack_edirect - unpack ediret executables

=head1 SYNOPSIS

  # Creating a new tool script
  use Rotifer::Tools::Rrsw;
  my $rrsw = Rotifer::Tools::Rrsw->create("unpack_edirect", { name => $target });
  $rrsw->get(@ARGV);

=head1 DESCRIPTION

Rotifer::Tools::Rrsw::unpack_edirect unpacks EDirect's executables.

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::Tools::Rrsw::unpack_edirect;

use strict;
use warnings;
use autodie qw(:all);
use Carp::Clan qr/^Rotifer::Tools::Rrsw/;
use File::Basename;
use File::Copy;
use File::Find;
use File::Which;
use IPC::Run qw(run);
use Moose;
use Rotifer;
extends 'Rotifer::Tools::Rrsw::unpack';

=head2 target

 Usage   : Rotifer::Tools::Rrsw::unpack_edirect->new(target => $target)
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

 Usage   : Rotifer::Tools::Rrsw::unpack_edirect->new(options => { gzip => "/usr/bin/gzip" })
 Function: module configuration
 Value   : hash reference as returned by Config::General
 Builder : _default_options

=cut

# Inherited from Rotifer::Tools::Rrsw::unpack

=head2 executable

 Usage   : Rotifer::Tools::Rrsw::unpack_edirect->new(executable => { gzip => "/usr/bin/gzip" })
 Function: path to tar, gzip and other executable
 Value   : hash reference with program names as keys
 Builder : _default_executable (uses File::Which)

=cut

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
    my $status = $self->next::method($source, $unpack);

    # Format unpacked databases
    mkdir($unpack) if ( ! -d $unpack );
    chdir($unpack);

    # Remove GO source files and non-executables
    if ($self->test) {
    } else {
      unlink("README");
      chdir("edirect");
      map { unlink $_ if ( -f $_ ) } ("Mozilla-CA.tar.gz","README","setup.sh","setup-deps.pl",glob("*.go"));
      chdir($unpack);
      map { move($_,basename($_)) } glob("edirect/*");
      chmod(0755, glob("*.Linux"));
      rmdir("edirect");
    }

    chdir($Rotifer::IWD);
    return $status;
}

__PACKAGE__->meta->make_immutable;
1;
