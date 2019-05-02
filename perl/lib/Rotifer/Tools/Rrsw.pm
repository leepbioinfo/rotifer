# POD documentation - main docs before the code

=head1 NAME

Rotifer::Tools::Rrsw - tools to update mirrors

=head1 SYNOPSIS

  # Creating a new parser
  use Rotifer::Tools::Rrsw;
  my $plugin = Rotifer::Tools::Rrsw->new("rsync", target => "blast");
  $plugin->execute($origin, $processed);

=head1 DESCRIPTION

Rotifer::Tools::Rrsw chooses and instantiates modules for
different download protocols and processing modules.

Each of these modules are controlled by configuration
files and should ownload and update a local copy of data.

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 DEPENDENCIES

=over

=item Carp::Clan
=item Module::Load

=back

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=head2 new

 Title   : new
 Usage   : Rotifer::Tools::Rrsw->new("arch",%newargs)
 Function: create a new instance of a plugin module
 Returns : whatever the module's new() method returns
 Args    : (string) module/plugin name and
	   (hash) arguments to the plugin's "new" method

	   For most Rrsw plugins, the "target" attribute
	   of the new() method is obligatory.

=cut

package Rotifer::Tools::Rrsw;

use Module::Load; # 'load' method
use Carp::Clan qr/^Rotifer::Parser/;
use Rotifer::Utils qw(describe_subclasses);
use strict;
use warnings;

# Set my subclass namespace
sub new {
    my ($self,$module,@args) = @_;
    confess "Unable to guess module!" if (!defined $module);
    if ($module eq "help") {
	print describe_subclasses(__PACKAGE__);
	exit 1;
    }
    $module = "${self}::$module";
    load "$module" || confess "Unable to load $module!";
    return $module->new(@args);
}

1;
