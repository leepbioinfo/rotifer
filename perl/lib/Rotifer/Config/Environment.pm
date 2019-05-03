=head1 NAME

Rotifer::Config::Environment: manages %ENV for Rotifer apps.

=head2 DESCRIPTION

This module provides configures Rotifer's environment variables.

=head3 Configuration

This module supports loading values from Rotifer's 
configuration file located at

 <rotifer install root>/etc/environmet.yml

or, if available,

 <rotifer install root>/etc/environmet.<hostname>.yml

where <hostname> is the name of the system were you are running
Rotifer without its domain name part, i.e. the output of the
command

 hostname | cut -f 1 -d "."

The format of the file is YAML (http://www.yaml.org). See examples
in this package for details.

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Rotifer::Config::Environment;

use Carp::Clan qr/^Rotifer::Config/;
use Config::Any;
use Cwd;
use Env qw(@PATH @LD_LIBRARY_PATH $ROTIFER_ROOT $ROTIFER_DATA);
use File::Basename;
use FindBin;
use lib "$FindBin::RealBin/../perl/lib";
use Moose;
use Sys::Hostname;
use strict;
use warnings;
use vars qw($__RealBasename $__RealBaseDir $__RealBaseEtc $__HostName);

=head2 SYSTEM ACCESSORS

Read-only methods to access your system information on
paths, environment variable, etc.

=cut

# Run as soon as the Perl compiler gets here!!!!
BEGIN {
    # Dirty trick to get the right paths when playing with 'perl -MRotifer -e'
    sub _get_basedir_from_lib_path {
	my $bdir = dirname($INC{"Rotifer.pm"});
	$bdir =~ s|/lib$||;
	$bdir =~ s|/perl$||;
	return Cwd::realpath($bdir);
    }

    # Basic settings
    $__RealBasename = basename($FindBin::RealScript);
    $__RealBaseDir  = $__RealBasename eq "-e" ? &_get_basedir_from_lib_path() : Cwd::realpath("$FindBin::RealBin/../");
    $__RealBaseEtc  = File::Spec->catfile($__RealBaseDir, qw(etc rotifer));
    ($__HostName = hostname) =~ s/\.\S+//;
};

=head2 host_base_name

 Title   : host_base_name
 Usage   : $env->host_base_name
 Function: name of the current host WITHOUT its domain name
 Returns : string
 Args    : none

=cut

sub host_base_name {
    return $__HostName;
}

=head2 real_basename

 Title   : real_basename
 Usage   : $env->real_basename
 Function: basename of your script, all symlinks resolved
 Returns : string
 Args    : none

=cut

sub real_basename {
    return $__RealBasename;
}

=head2 real_basedir

 Title   : real_basedir
 Usage   : $env->real_basedir
 Function: actual path to the root of Rotifer's
           installation (all links resolved) 
 Returns : string
 Args    : none

=cut

sub real_basedir {
    return $__RealBaseDir;
}

=head2 real_system_config_root

 Title   : real_system_config_root
 Usage   : $env->real_config_root
 Function: actual path to the root of Rotifer's
           configuration directory tree
 Returns : string
 Args    : none

=cut

sub real_system_config_root {
    return $__RealBaseEtc;
}

=head2 real_system_data_root

 Title   : real_system_data_root
 Usage   : $env->real_system_data_root("/usr/loca/data")
 Function: get/set path to Rotifer's data directories
 Returns : list of directories
 Args    : none
 Default : <rotifer installation directory>/data

=cut

sub real_system_data_root {
    my $self = shift if (defined $_[0] && ($_[0] eq __PACKAGE__ || ref($_[0]) eq __PACKAGE__));
    $ROTIFER_DATA = shift if (scalar @_);
    return $ROTIFER_DATA;
}

=head2 MODIFIERS

=head2 add_to_path

 Title   : add_to_path
 Usage   : $env->add_to_path("~/bin")
 Function: expand the executables search path
 Returns : string (new concatenated list)
 Args    : list of directories to be added

=cut

sub add_to_path {
    my $self = shift if (defined $_[0] && ($_[0] eq __PACKAGE__ || ref($_[0]) eq __PACKAGE__));
    my @dir = @_;
    if (scalar @dir) {
	my %path = map { ($_,1) } @PATH;
	push(@PATH, grep { !exists $path{$_} } @dir);
    }
    return $ENV{"PATH"};
}

=head2 change_environment

 Title   : change_environment
 Usage   : $env->change_environment(\%o, A => 1, PATH => ["~/bin"])
 Function: set values of environment variables
 Args    : hash with new values (key: variable name)
 Returns : nothing
 Options : hash reference (\%o). Available options are

  name    : description                     : type    : default
 ---------:---------------------------------:---------:---------
  create  : create new variables            : boolean : 1 
  push    : for lists, use Perl's "push"    : boolean : 0 
  unshift : for lists, use Perl's "unshift" : boolean : 0 

 These options apply only to variables that can hold list of
 values (e.g. PATH and LD_LIBRARY_PATH). Such variables are
 set using references to arrays.

=cut

{
    no strict "refs";

    sub change_environment {
	my $__self = shift if (defined $_[0] && ($_[0] eq __PACKAGE__ || ref($_[0]) eq __PACKAGE__));
	my $__opts = { %{ shift @_ } } if (defined $_[0] && UNIVERSAL::isa($_[0],"HASH"));
	$__opts->{create} = 1 unless (exists $__opts->{create});
	my %__hash = @_;

	# Process each variable
	foreach my $__name (keys %__hash) {
	    next unless ($__opts->{create} || $ENV{$__name});
	    my $__value = $__hash{$__name};

	    if (ref $__value eq "ARRAY") {
		Env->import("\@$__name");
		if ($__opts->{push}) {
		    push(@{"$__name"},@$__value);
		} elsif ($__opts->{unshift}) {
		    unshift(@{"$__name"},@$__value);
		} else {
		    @{"$__name"} = @$__value;
		}
	    }

	    else {
		Env->import("\$$__name");
		$$__name = $__value;
	    }
	}

	return ;
    }

    use strict 'refs';
};

=head2 PRIVATE METHODS

=head2 _load_config

 Title   : _load_config
 Usage   : $env->_load_config()
 Function: parse and load <rotifer_root>/etc/config/environment.yml
 Returns : string
 Args    : none

=cut

sub _load_config {
    # Laod configuration file
    my $BaseConfig = {};
    my $host = &host_base_name();
    my $file = "$__RealBaseEtc/config/environment.yml";
    $file = "$__RealBaseEtc/config/environment.${host}.yml" if ( -e "$__RealBaseEtc/config/environment.${host}.yml" );
    return unless ( -r $file );
    $BaseConfig = Config::Any->load_files({ files => [ $file ],
					    use_ext => 1,
					    flatten_to_hash => 1,
					  });
    change_environment({ push => 1, create => 1 }, %{ $BaseConfig->{$file} });
};

# Run as soon as the Perl compiler gets here!!!!
BEGIN {
    # Parse configured environment
    _load_config();

    # Make sure basics are always there
    $ROTIFER_ROOT = $__RealBaseDir unless (defined $ROTIFER_ROOT);
    real_system_data_root("$__RealBaseDir/data") if (!exists $ENV{ROTIFER_DATA});

    # Add suport for user's settings, data, etc
    add_to_path("$ENV{HOME}/.rotifer/bin");
};

# Make perl compiler happy returning meaningless true value
__PACKAGE__->meta->make_immutable;
1;
