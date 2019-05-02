=head1 NAME

ROTIFER - Rapid Open source Tools and Infrastructure For
          data Exploration and Research

=head2 DESCRIPTION

Rotifer is a somewhat arbitrary collection of tools for data analysis.
It is also an infrastructure for the development of new data analysis 
and data management tools. Technically, Rotifer intends to provide a 
high-level application development framework and a collection of tools
based on this framework.

In pratical terms it is just a collection of programs, scripts and 
libraries I wrote for my bioinformatics research and latter compiled
in a somewhat user-friendly package. 

=head2 FEATURES

=over

=item * A set of programs for processing and analysing data in tabular format

=item * Standard, simplified interfaces for commonly used command line bioinformatics tools

=item * An infrastructure for users to manage and expand their Rotifer installation

=item * Libraries/modules to easy development of new Rotifer tools 

=back

=head2 COMPONENTS

Rotifer contains several scripts organized hierarchically in sections and by the
assignment of tags. Users may add new scripts using Rotifer's mantras system.

=head3 The rotifer command

The entry point to get to know Rotifer is the rotifer command. It provides access
to all the documentation on the systems including comments, examples and tags added
by the user. It also provides controls for the installation/removal of system 
components and its configuration. You may access rotifer's documentation online or
by typing 

 rotifer help

if Rotifer is already installed in your system.

=head3 RotIFeR libraries

 Rotifer
 Rotifer::Config::Environment

=head2 GOALS

Rotifer's utmost goal is to be B<useful>. If this objective is not achieved then I will just
have wasted my time writing it. Another goal is B<simplicity>, i.e. B<easy> of use. The latter
goal often contradicts the third goal which is to B<avoiding reimplementation> of the same 
solution but I try to devote the same amount of energy to both problems. Unfortunately, 
it seems that only the first goal is easy to verify by just checking whether there are any
results. A program is only simple after a user learns how to use it and every user should
understand that Rotifer has its own culture. I hope to have made learning this culture both
effortless and enjoyable.

The implementation of Rotifer tries to avoid is the (re)creation of dozens of tools that perform 
essentially the same task but with minor variations. In Rotifer I almost always choose to freeze
each solution in a B<flexible> and stable tool and use Rotifer's own capacity to store and 
document the various ways a tool is used to solve similar problems without requiring the tool 
to be rewritten.

=over

=item * Make it easy to create programs with the same command line behaviour, i.e. 
(1) programs that parse the command line in a semi standardized way and
(2) that provide the same 

=back

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Rotifer;

use Cwd qw(abs_path getcwd);
use Config::Any;
use File::Basename;
use FindBin;
use lib "$FindBin::RealBin/../perl/lib";
use Rotifer::Config::Environment;
use strict;
use warnings;

# Set values to some global variables
use vars qw($BaseConfig $RealBasename $RealBaseDir $RealScriptEtc);
our $VERSION = "0.86";
our ($IWD) = abs_path(getcwd());

# Default Rotifer configuration
BEGIN {
    # Base paths
    $RealBasename  = Rotifer::Config::Environment->real_basename;
    $RealBaseDir   = Rotifer::Config::Environment->real_basedir;
    my $etc = Rotifer::Config::Environment->real_system_config_root;
    $RealScriptEtc = "$etc/$RealBasename";
}

=head2 PATH ACCESSORS

General methods that encapsulate where Rotifer components search for
directories, files, etc. 

=head2 sysdir

 Title   : sysdir
 Usage   : $env->sysdir("data")
 Function: get path to directories based on your Rotifer installation
 Returns : path to a file or directory
 Args    : (string) optional subdirectory path

=cut

sub sysdir {
    my $self = shift if (defined $_[0] && ($_[0] eq __PACKAGE__ || ref($_[0]) eq __PACKAGE__));
    my $path = File::Spec->catfile(Rotifer::Config::Environment->real_basedir);
    if (my @path = @_) {
	my $subdir = shift(@path);
	unshift(@path,"rotifer") if ((!defined $path[0] || $path[0] ne "rotifer") && grep { $subdir eq $_ } qw/etc lib share/);
	$path = File::Spec->catfile($path, $subdir, @path);
    }
    return $path;
}

=head2 userdir

 Title   : userdir
 Usage   : $env->userdir("data")
 Function: get path to user-specific Rotifer directories
 Returns : path to a file or directory
 Args    : (string) optional subdirectory path

=cut

sub userdir {
    my $self = shift if (defined $_[0] && ($_[0] eq __PACKAGE__ || ref($_[0]) eq __PACKAGE__));
    my $path = File::Spec->catfile($ENV{HOME},".rotifer");
    $path = File::Spec->catfile($path, @_) if (scalar @_);
    return $path;
}

=head2 config

 Title   : config
 Usage   : Rotifer->config
 Function: contents of etc/rotifer.yml
 Returns : hash reference
 Args    : none

=cut

sub config {
    if (!defined $BaseConfig) {
	# Choose global configuration file
	my $config = Rotifer::Config::Environment->host_base_name;
	my $etc = &sysdir("etc");
	$config = "$etc/rotifer.${config}.yml" if (defined $config);
	$config = "$etc/rotifer.yml" if ( ! -e $config );

	# Parse configuration file
	$BaseConfig = -e $config
	    ? Config::Any->load_files({ files => [ $config ],
					use_ext => 1,
					flatten_to_hash => 1,
				      })->{$config}
	: {};
    }

    return $BaseConfig;
}

# Make perl compiler happy returning meaningless true value
1;
