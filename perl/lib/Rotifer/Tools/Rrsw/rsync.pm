# POD documentation - main docs before the code

=head1 NAME

Rotifer::Tools::Rrsw::rsync - download/update mirror using rsync 

=head1 SYNOPSIS

  # Creating a new tool script
  use Rotifer::Tools::Rrsw;
  my $rrsw = Rotifer::Tools::Rrsw->create("rsync", { name => $target });
  $rrsw->execute($origin, $processed);

=head1 DESCRIPTION

Rotifer::Tools::Rrsw::rsync is a wrapper for the rsync program
that follows the interface required to be a component of rrsw.

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::Tools::Rrsw::rsync;

use strict;
use warnings;
use autodie qw(:all);
use Carp::Clan qr/^Rotifer::Tools/;
use File::Which;
use IPC::Run qw(run);
use Moose;
use Rotifer;
use Scalar::Util qw(blessed);

=head2 debug

 Usage   : Rotifer::Tools::Rrsw->new(debug => 1)
 Function: activate debbuging messages
 Value   : string
 Required: 1

=cut

has 'debug' => (
    is       => "rw", 
    isa      => "Int",
    default  => 0,
    );

=head2 target

 Usage   : Rotifer::Tools::Rrsw->new(target => $target)
 Function: basename of the target
 Value   : string
 Required: 1

=cut

has 'target' => (
    is       => "rw", 
    isa      => "Str",
    required => 1,
    );

=head2 executable

 Usage   : Rotifer::Tools::Rrsw->new(executable => "/usr/bin/rsync")
 Function: path to rsync's executable
 Value   : string
 Builder : _default_executable (uses File::Which)

=cut

has 'executable' => (
    is       => "rw", 
    isa      => "Str",
    lazy     => 1,
    builder  => '_default_executable',
    );

sub _default_executable {
    my $rsync = which("rsync");
    die "Could not find the rsync executable" unless (defined $rsync);
    return $rsync;
}

=head2 options

 Usage   : Rotifer::Tools::Rrsw->new(options => { "rsync" => [ qw/--partial/ ] })
 Function: command line options for rsync
 Value   : array reference
 Builder : _default_options
 Trigger : _options_trigger
=cut

has 'options' => (
    is       => "rw", 
    isa      => "HashRef[ArrayRef]",
    lazy     => 1,
    builder  => '_default_options',
    trigger  => \&_options_trigger,
    );

sub _default_options {
    return { "rsync" => [ "-arHv", "--delete", "--prune-empty-dirs", "--delete-excluded" ] };
}

sub _options_trigger {
    my ($self, $new, $old) = @_;
    if (defined $new && !scalar(keys %$new)) {
	$self->options($self->_default_options);
    }
}

=head2 test

 Usage   : Rotifer::Tools::Rrsw->new(test => 1)
 Function: run in test mode
 Value   : Boolean

=cut

has 'test' => (
    is       => "rw", 
    isa      => "Bool",
    lazy     => 1,
    default  => 0,
    );

=head2 METHODS

=head2 execute

 Title   : execute
 Usage   : $count = $rrsw->execute($orig, $dest)
 Function: download/sync mirror
 Returns : list of files donwloaded
 Args    : origin and destination

=cut

sub execute {
    my ($self, $orig, $dest) = @_;

    # Prepare command-line
    my $name = $self->target;
    my @cmd = ($self->executable, @{$self->options->{rsync}});
    push(@cmd, $self->test ? "-n" : "--no-motd");
    push(@cmd, "-f", ". $Rotifer::RealScriptEtc/${name}.rules") if ( -f "$Rotifer::RealScriptEtc/$name.rules" );
    push(@cmd, $orig, $dest);

    # Run!
    print STDERR join(" ",@cmd),"\n" if ($self->debug > 3);
    run(@cmd);
    return !$?;
}

__PACKAGE__->meta->make_immutable;
1;
