# POD documentation - main docs before the code

=head1 NAME

AppConfig::AutoDoc - parse options from command
line, configuration files and/or CGI parameters

=head1 SYNOPSIS

  # Minimal (not very useful)
  use AppConfig::AutoDoc;
  my $config = AppConfig::AutoDoc->new();
  $config->parse_cmdline_files();

  # Creating a boolean flag and using it
  #
  use AppConfig::AutoDoc qw(:argcount);
  my $config = AppConfig::AutoDoc->new();
  $config->define("remove_old" => {
                  ALIAS    => "r",
                  ARGCOUNT => ARGCOUNT_NONE,
                  SUMMARY  => "remove old fastas",
                  DEFAULT  => 0,
                 });
  $config->parse_cmdline_files();
  unlink(<*.fa>) if ($config->remove_old);

=head1 DESCRIPTION

This module allows easy definition and standadized parsing of program
options from the command line, CGIs and/or configuration file(s). It
derives most of its functionality from L<AppConfig> module by Andy
Wardley.

It extends AppConfig by adding some predefined options and
auto-generated usage messages and documentation. It also extends
AppConfig's functionality by introducing new tags to the L<defline>
and L<new> methods (see below).

See "perldoc L<AppConfig>" for details.

=head1 Forbidden option names

Since AppConfig and Appconfig::AutoDoc use AUTOLOAD to provide
accessors to configured options, any option whose name conflicts
with the name of a subroutine in AppConfig or AppConfig::AutoDoc
can't be differentiated from the subroutine and will result in the
program aborting execution and dumping an error message.

In order to avoid such conflicts, most methods in this library are
named starting with an underscore ("_") but a few subroutines names
are not prefixed with _ and may cause conflict. The non-dashed
subroutines currently implemented in AppConfig::AutoDoc are:

=over

=item new

=item get_options

=item get_type

=item get_aliases

=item get_defaults

=item define

=item to_file

=item to_string

=item option_to_string

=item parse_cmdline_files

=item getargs

=back

The user must also avoid method names that might conflict with
method names from the following libraries

=over

=item AppConfig

=item File::Basename

=item Getopt::Long

=item Exporter

=item UNIVERSAL

=back

=head1 Automatically generated options

AppConfig::AutoDoc automatically defines the following options

 --help or --usage or -h
 --debug
 --configfile
 --configdump

Options --usage triggers execution of internal routines that
automatically generate and print documentation added by the user
to its program in POD format and derived by AppConfig::AutoDoc
from the program option's definitions.

Option --configfile directs the parser to process one or more
configuration files in AppConfig's .INI style format (see perldoc
AppConfig). Parsing of these files is automatic and more than one
configuration file may be chosen by repeating this option in the
command line (--configfile file1 --configfile file2). Note that
this option will not be available if there are no user defined
options.

Option --debug does not trigger execution of any predefined method
unless its argument is 'config' (--debug config). In this case,
AppConfig::AutoDoc formats all information about its internal
state and configuration parameter values and prints it to STDERR.

For more details about predefined options, see the documentation of
each triggered method in the APPENDIX section below.

=head1 AUTHOR - Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package AppConfig::AutoDoc;

use strict;
use warnings;

use Carp::Clan qr/^AppConfig::AutoDoc/;
use File::Basename;
use Getopt::Long;
use Exporter;
use AppConfig qw(:expand :argcount);
use base qw(AppConfig Exporter);
use vars qw($AUTOLOAD %DEFAULT_OPTIONS);

our $VERSION = "0.1";

# Expanding AppConfig's exported symbols list
our @EXPORT_OK   = (@AppConfig::EXPORT_OK);
our %EXPORT_TAGS = (%AppConfig::EXPORT_TAGS);

# The hash of hashes below defines the default options
# that are available from this class. For details, see
# 'perldoc AppConfig' and the 'define' method below
our %DEFAULT_OPTIONS = (# Configuration file
    'configfile' => {
        ARGCOUNT => ARGCOUNT_LIST,
        DEFAULT  => [],
        SUMMARY  => "configuration file(s) to parse",
    },

    'configdump' => {
        ARGCOUNT => ARGCOUNT_NONE,
        DEFAULT  => undef,
        SUMMARY  => "Print active options to standard output as configuration file.",
    },

    # Debug level or procedure
    'debug' => {
        DEFAULT  => 0,
        ARGCOUNT => ARGCOUNT_ONE,
        SUMMARY  => "set debug level (integer or 'config')",
    },

    # Short command line summary
    'help' => {
        ALIAS    => 'usage|h',
        DEFAULT  => 0,
        ARGCOUNT => ARGCOUNT_NONE,
        SUMMARY  => "show program documentation"
    },

    # Print version
    'version' => {
        DEFAULT  => 0,
        ARGCOUNT => ARGCOUNT_NONE,
	ACTION   => sub {
		use Rotifer;
		use Rotifer::Utils qw(pod2summary);
		use File::Basename;
		my @version = ("$0:");
		my @desc = pod2summary({ section => ["AUTHOR"] }, $0);
		push(@version,$desc[0][1]);
		my $versionNumber = defined($main::VERSION) ? $main::VERSION : $Rotifer::VERSION;
		push(@version,"Version: $versionNumber");
		push(@version,"Authors: $desc[0][2]");
		print join("\n",@version),"\n";
		exit 1;
	},
        SUMMARY  => "show program version"
    },
);

# Default configuration parameters for new()
our %DEFAULT_PARAMS = (# Default error handler
    'ERROR' => \&_default_error_handler,
);

# Extended attributes for new() and defined()
our %EXTENDED_CLASS_TAGS = (# Expected behaviour for @ARGV
    'DEACTIVATE_EXTENSION' => [],
    'ERASE_DEFAULTS'       => 1,
    'EXPECT_INPUT'         => 1,
    'EXPAND_ARGV'          => 0,
    'EXPAND_STDIN'         => 0,
    'IGNORE_UNKNOWN'       => 0,
);
our %EXTENDED_OPTION_ATTRIBUTES = (# Short option description
    'SUMMARY' => "",
    'SECTION' => "",
);

# Controlling automatic file parsing for arrays, hashes, etc.
our %EXPAND_FILES_OPTIONS = (# How to parse
    'COLUMN_DELIMITER' => '\s+',
);

=head2 Object constructor

=head2 new

 Title   : new
 Usage   : $config = $self->new()
 Function: create object
 Returns : nothing
 Args    : hash reference (optional) and/or
           hash of hashes defining options (optional)

 Note: AppConfig allows setting some options by passing
       a hash reference as the first argument to L<new>
       (preceding any definition for options). This
       module supports the same tags as AppConfig and
       extends it by adding the following set of tags:

       EXPECT_INPUT => 0 or 1 (boolean, default is 1)

       The EXPECT_INPUT tag controls whether or not program
       execution should be ABORTED IF the program receives
       NO ARGUMENTS, i.e. if @ARGV is empty program options
       are parsed AND the program is NOT receiving input from
       a pipe. This feature is enabled by default.

       EXPAND_ARGV => 0 or 1 (boolean, default is 0)

       The EXPAND_ARGV tag is a boolean and, therefore, can take
       values 0 (false) or 1 (true). If set to true (1) this tag
       will activate replacement of all file names listed in
       @ARGV by the contents of each row of each file. Such
       processing will be performed only after parsing (and
       removal) of all command line options and parameters from
       @ARGV and permits adding to @ARGV lists of arguments much
       longer than the limits set by the shell.

       EXPAND_STDIN => 0 or 1 (boolean, default is 0)

       The EXPAND_STDIN tag is also boolean and analogous to
       EXPAND_ARGV but, instead of placing the contents of a file
       named in @ARGV, it transfers to @ARGV the rows in any
       text message passed to the program's standard input through
       pipe's or other means.

       IGNORE_UNKNOWN => 0 or 1 (boolean, default is 0)

       By default, AppConfig::AutoDoc aborts parsing of @ARGV and
       prints a detailed error message whenever it encounters a
       command line option that was not defined by the program's
       author. Command line options are identified by starting with
       '-' or '--' followed by an arbitrary string. Setting the flag
       IGNORE_UNKNOWN to 1 (true) makes the parser bypass unknown
       options, the same way it does with other arguments that
       remain in @ARGV after parsing of the command line.

       ERASE_DEFAULTS => 0 or 1 (boolean, default is 1)

       The ERASE_DEFAULTS flag controls whether DEFAULT values set to
       hashes, arrays and the --configfile option in the application
       should be erased before parsing of the command line and
       configuration files. Setting this flag to FALSE makes all
       new values set the application's user to be appended to hashes
       and arrays in the configuration object, instead of entirely
       overwriting these lists.

       DEACTIVATE_EXTENSION => reference to an array (default: empty)

       By default, when at least one command line option is defined,
       AppConfig::AutoDoc automatically adds its pre-defined options
       (help, debug and configfile). By giving a list of option names
       as elements of an anonymous array to this parameter, the user
       may disable any of these pre-defined options.

       =-----------------------------------------------------=

       In addition to the hash reference above, one can also modify
       the definition of pre-configured options introduced by this
       class. To achieve this, just add arguents to define these
       options as if they were not pre-defined one, but ommit
       whatever parameters you do not want to change. For example,

       "configfile" => { DEFAULT => [ "/etc/app.conf" ] }

       would make the file /etc/app.conf the default configuration
       file for your program.

       See the section on the method "define()" for details on how to
       set up options.

=cut

sub new {
    my ($class, @args) = @_;

    # Copy default values for configuration hash
    my $opt = { %DEFAULT_PARAMS };
    my %added_tags = %EXTENDED_CLASS_TAGS;

    # Copy control parameters but isolate extended
    # tags to avoid breaking AppConfig::new()
    if (ref($args[0]) eq "HASH") {
        my $parameters = shift(@args);
        foreach my $param (keys %{$parameters}) {
            if (exists $EXTENDED_CLASS_TAGS{$param}) {
                $added_tags{$param} = $parameters->{$param};
            } else {
                $opt->{$param} = $parameters->{$param};
            }
        }
    }

    # Call method new from parent class (AppConfig)
    my $self = $class->SUPER::new($opt);
    map { $self->{$_} = $added_tags{$_} } keys %added_tags; # Store extended tags

    # Copy and change default parameters to match user input
    my %args    = @args;
    my %options = %DEFAULT_OPTIONS;
    if (scalar @args) {
        foreach my $option (keys %options) {
            next unless (exists $args{$option});
            foreach my $param (keys %{ $args{$option} }) {
                next if ($param eq 'ARGCOUNT'); # Users cannot change the type of predefined options
                $options{$option}->{$param} = $args{$option}->{$param};
            }
            delete $args{$option};
        }
    }

    # Create default options (remove configfile if there are no user defined options)
    delete $options{'configfile'} if (!scalar(keys %args));
    foreach my $option (@{ $added_tags{'DEACTIVATE_EXTENSION'} }) {
        delete $options{$option};
    }

    # Create options
    $self->define(%options);

    # Create user requested options
    $self->define(%args);

    # Create AppConfig's object and define our default options
    return $self;
}

=head2 Accessors designed for this class

=head2 get_options

 Title   : get_options
 Usage   : $self->get_options()
 Function: get the list of all options for this program
 Returns : array of strings
 Args    : none

=cut

sub get_options {
    my $self = shift;
    my %hash = $self->varlist(".");
    return sort keys %hash;
}

=head2 get_type

 Title   : get_type
 Usage   : $config->get_type()
 Function: get the ARGCOUNT attribute of a variable
 Returns : string
 Args    : string

=cut

sub get_type {
    my ($self, $var) = @_;

    # Return aliases to one variable
    return undef if (!defined $var);
    my $type = undef;
    eval { $type = $self->_argcount($var) };
    die "Problem accessing definition for option $var:\n$@" if ($@);
    TYPE: {
        (!defined $type) && confess "Undefined variable: $var";
        ($type == 0) && do { $type = "boolean"; last TYPE; };
        ($type == 1) && do { $type = "scalar";  last TYPE; };
        ($type == 2) && do { $type = "array";   last TYPE; };
        ($type == 3) && do { $type = "hash";    last TYPE; };
    };
    return $type;
}

=head2 get_aliases

 Title   : get_aliases
 Usage   : $config->get_aliases()
 Function: get a list of alias for options
 Returns : list (if an option name is provided, might be empty) or
           hash (option name => reference to array of get_aliases)
 Args    : string

=cut

sub get_aliases {
    my ($self, $var) = @_;

    # Note: Unfortunately, AppConfig does not provide an
    #       accessor for its aliases definitions, therefore
    #       I'm forced to make direct access to its internal
    #       data structure: if that changes without warning,
    #       this method will fail
    my $aliases = $self->_state->{ALIASES};

    # Return aliases to one variable
    if (defined $var) {
        return () if (!exists $aliases->{$var});
        return @{$aliases->{$var}};
    }

    # Return hash of all aliases
    my %hash = ();
    foreach my $varname (keys %$aliases) {
        $hash{$varname} = [ @{$aliases->{$varname}} ];
    }
    return %hash;
}

=head2 get_defaults

 Title   : get_defaults
 Usage   : $config->get_defaults("help")
 Function: get default values for a variable
 Returns : depends on the ARGCOUNT flag used when defining the option
           ARGCOUNT_NONE => boolean (scalar set to either 0 or 1)
           ARGCOUNT_ONE  => scalar (string)
           ARGCOUNT_HASH => reference to hash
           ARGCOUNT_LIST => reference to array
           undef         => when no default values were set
 Args    : string, option name

=cut

sub get_defaults {
    my ($self, $var) = @_;
    # Note: Unfortunately, AppConfig does not provide an
    #       accessor for its aliases definitions, therefore
    #       I'm forced to make direct access to its internal
    #       data structure: if that changes without warning,
    #       this method will fail
    my $defaults = $self->_state->{DEFAULT};
    return undef if (!defined $var);
    return undef if (!exists $defaults->{$var});
    $defaults = $defaults->{$var};
    return $defaults;
}

=head2 Methods that override AppConfig's methods

=head2 define

 Title   : define
 Usage   : $self->define("size" => { 'ALIAS' => 's', DEFAULT => 1 })
 Function: add options to your program's interface
 Returns : true or false
 Args    : hash of hashes

 Note: This method extends AppConfig's define method to provide
       documentation capabilities. See L<AppConfig::define>
       for base functionality and use the following extra tags
       added by this modules:

       SUMMARY => very short description of the option  (string)
       SECTION => string used to group related options in the
                  auto-generated POD documentation

       In addition to the new tags above, this module also
       introduces a set of pre-defined ACTIONs that may be
       executed every time an option is set. To activate a
       pre-defined action, set the value of the ACTION key to
       the name of the desired action. Example

       $config->define('option' => { ACTION => 'EXPAND_FILES' })

       List of pre-defined actions:

       'NO_EXPECT_INPUT' => turn off input monitoring when a
                            variable is set.

        Notice that this flag will evaluate the value set as either
        true or false so it is mostly useful for boolean options.

       'EXPAND_FILES' => parse and load the contents of files

        EXPAND_FILES effect will depend on the variable's
        ARGCOUNT value. If the variable type is ARGCOUNT_HASH,
        EXPAND_FILES will extract ONLY the first two columns of
        text files identified by the key "file". and treat the first column as keys and the
        second column as values for the internal hash. NOTE:
        when loading hashes, repeated values in the first
        column will keep ONLY THE LAST ROW VALUES. For
        ARGCOUNT_LIST, each row in a text file is loaded to a
        separate element of the option's array and for
        ARGCOUNT_ONE the entire file is loaded as a single
        string. This parameter has no effect if ARGCOUNT is
        set to ARGCOUNT_NONE.

        When requesting hash expansion from files, only the values
        are used, i.e. file names should be set like in this
        example:

        my_prog.pl --hash file=file1.txt --hash file=file2.txt

=cut

sub define {
    my ($self, %args) = @_;

    # Process extra flags in input hash
    my @expand = ();
    foreach my $option (keys %args) {
        if (__PACKAGE__->can($option)) {
            my $message = "$option is a true subroutine in AppConfig::AutoDoc!!!

            You can't use option names that conflict with internal method names (see 'perldoc AppConfig::AutoDoc')";
            $self->_default_error_handler($message);
        }

        # Managing extended attributes
        foreach my $tag (keys %EXTENDED_OPTION_ATTRIBUTES) {
            $self->_extended_attribute($option, $tag, $args{$option}->{$tag});
            delete $args{$option}->{$tag};
        }

        # Dealing with pre-defined ACTIONs
        if (exists $args{$option}->{'ACTION'} && !ref $args{$option}->{'ACTION'}) {
            foreach my $action (split(/\,+/,$args{$option}->{'ACTION'})) {
                if ($action eq 'EXPAND_FILES') {
                    $args{$option}->{'ACTION'} = \&_expand_file;
                    push(@expand, $option);
                }
                elsif ($action eq 'NO_EXPECT_INPUT') {
                    $args{$option}->{'ACTION'} = \&_no_expect_input;
                }
            }
            die "Proposed ACTION for option $option won't compile: ".$args{$option}->{'ACTION'}." "
            if (exists $args{$option}->{'ACTION'} && !ref $args{$option}->{'ACTION'});
        }
    }

    # Executing the original define subroutine
    my $ret = $self->SUPER::define(%args);

    # AppConfig has yet another bug: it does not apply ACTION to default values
    # Let's make yet another ugly hack to fix this for our predefined actions
    foreach my $option (@expand) {
        if ($self->_argcount($option) == ARGCOUNT_HASH) {
            &_expand_file($self->_state, $option);
        } elsif ($self->_argcount($option) == ARGCOUNT_LIST) {
            &_expand_array($self->get($option));
        } else { # Scalar or boolean
            &_expand_file($self->_state, $option, $self->get($option));
        }
    }

    return $ret;
}

=head2 Methods that extend AppConfig's interface

=head2 to_file

 Title   : to_file
 Usage   : $config->to_file("file.conf")
 Function: dump configuration to a file
 Returns :
 Args    : (string) output file name
           (optional) list of options

 The generated file maybe parse by AppConfig::AutoDoc's
 "file" method to create a configuration object identical
 to the current one. By default, all options are included
 in the output.

=cut

sub to_file {
    my ($self, $file, @options) = @_;
    open(TOF,">$file") || die "Could not dump configuration to file $file";
    print TOF $self->to_string(@options);
    close(TOF);
    return 1;
}

=head2 to_string

 Title   : to_string
 Usage   : $config = $config->to_string()
 Function: dump all non-default options as a string
 Returns : string
 Args    : (optional) list of options

 The returned string maybe parse by AppConfig::AutoDoc's
 "file" method to set a configuration object to identical
 option values. By default, all options are included in the
 output.

=cut

sub to_string {
    my ($self, @options) = @_;
    @options = $self->get_options if (!scalar @options);

    # Header
    my $text = "# Configuration file generated by $0 (PID: $$)\n\n";

    # Identify sections and isolate them
    my @main     = (); # Main options
    my $sections = {}; # Options from sections
    foreach my $option (@options) {
        next if ($option eq 'configfile'); # This one makes no sense to dump since only getargs uses it
        my $type = $self->get_type($option);
        if ($option =~ /_/ && $type ne 'boolean') { # Might belong to a section
            my ($section,$name) = split(/_/,$option,2);
            push(@{ $sections->{$section} }, $option);
        } else {
            push(@main, $option);
        }
    }

    # Append main section
    $text .= "# Main section\n";
    foreach my $option (@main) {
        $text .= $self->option_to_string($option);
    }

    # Append other sections
    # Lets not start sections until there is more than one option to print
    $text .= "\n# Other sections\n";
    my @names = sort { scalar(@{$sections->{$a}}) <=> scalar(@{$sections->{$b}}) || $a cmp $b } keys %$sections;
    foreach my $secname (@names) {
        my @section = ("\n[$secname]\n");
        foreach my $name (sort @{ $sections->{$secname} }) {
            my $string = $self->option_to_string($name, 1);
            push(@section, $string) if (defined $string && length($string));
        }
        $text .= join('',@section) if (scalar(@section) > 1);
    }

    return $text;
}

=head2 option_to_string

 Title   : option_to_string
 Usage   : $string = $config->option_to_string('debug')
 Function: retrun a string representing an option
 Returns : string
 Args    : (string)  option name
           (integer) omit section name

 The returned string maybe parse by AppConfig::AutoDoc's
 "file" method to reset the dumped option.

=cut

sub option_to_string {
    my ($self, $opt, $in_other_section) = @_;
    my $name = $opt;
    $name =~ s/^[^_]+_// if ($in_other_section);

    return undef unless ($self->_exists($opt));
    my $text = "";

    # Boolean
    if ($self->_argcount($opt) == 0) {
        $text .= $self->get($opt) ? "$name = 1\n" : "$name = 0\n";
    }

    # Scalar
    elsif ($self->_argcount($opt) == 1) {
        my $value = $self->get($opt);
        return "# $name = # Undefined values aren't supported by AppConfig::File!\n" if (!defined $value);
        return "# $name = # Zero length strings aren't supported by AppConfig::File!\n" if (!length $value);
        $value =~ s/\n/\\\n/g;
        $value = "\"$value\"" if ($value =~ /[\s\#]+/);
        $text .= "$name = $value\n";
    }

    # Array
    elsif ($self->_argcount($opt) == 2) {
        my $ref = $self->get($opt);
        foreach my $value (@$ref) {
            if (!defined $value) {
                $text .= "# $name = # Undefined values aren't supported by AppConfig::File!\n";
            } elsif (!length $value) {
                $text .= "# $name = # Zero length strings are not supported by AppConfig::File!\n";
            } else {
                $value = "\"$value\"" if ($value =~ /[\s\#]+/);
                $text .= "$name = $value\n";
            }
        }
    }

    # Hash
    elsif ($self->_argcount($opt) == 3) {
        my $ref = $self->get($opt);
        foreach my $key (keys %$ref) {
            my $value = $ref->{$key};
            if (!defined $value) {
                $text .= "$name $key = # Note: use of undefined for hash values is supported by AppConfig::File\n";
            } elsif (!length $value) {
                $text .= "$name $key = # Note: zero length strings for hash values are supported by AppConfig::File\n";
            } else {
                $value = "\"$value\"" if ($value =~ /[\s\#]+/);
                $text .= "$name $key = $value\n";
            }
        }
    }

    return $text;
}


=head2 parse_cmdline_files

 Title   : parse_cmdline_files
 Usage   : $config->parse_cmdline_files
 Function: parse command line arguments and
           configuration file(s) set by --configfile
 Returns :
 Args    : (optional) reference to an array to parse
                      instead of the default @ARGV

           (optional) string to generate help message
                      if this argument is given, the
                      program will dump the help
                      message to STDOUT and won't exit

 This method implicates the following options are added
 to your application interface:

 Option    : --configfile
 Aliases   :
 Arguments : a list of file name(s)
 Default   : off

=cut

sub parse_cmdline_files {
    my ($self, $args, @message) = @_;

    # Since AppConfig does not erase default entries while loading arrays
    # and hashes, we have to empty these lists here!!!! This is an AppConfig bug!!!!!
    #print STDERR join(" ","BEGINING:",map { my $v = $self->get($_); $v = join(",",(ref($v) eq 'HASH' ? %$v : ref($v) ? @$v : defined $v ? $v : '')); $_."=".$v } $self->get_options),"\n";
    $self->_config2hash($self->{"ERASE_DEFAULTS"});

    # Parsing command line arguments and (maybe) configuration file(s)
    $self->getargs(defined $args && scalar @$args ? $args : undef);
    my %cmdline = $self->_config2hash(); # Store AppConfig state
    my @options = $self->get_options;
    if (grep { $_ eq 'configfile' } @options) {
        if (!scalar @{ $self->get("configfile") }) { # Use default configuration files unless other were chosen by the user 
            my $default_configfiles = $self->get_defaults("configfile");
            my @default_configfiles = grep { -f $_ } @$default_configfiles; # Use default configuration files only if they exist!
            map { $self->set("configfile",$_) } @default_configfiles;
        }
        delete $cmdline{"configfile"} if (exists $cmdline{"configfile"}); #
        if (scalar @{$self->configfile}) {
            # Parsing configuration files (--configfile option)
            foreach my $file (@{ $self->configfile }) {
                die _stack_trace_dump("Error parsing command line: could not read configuration file $file") if (! -s $file);
                $self->file($file);
            }
            $self->_reset_from_hash(%cmdline); # Command line has precendence: enforce
            @options = $self->get_options;
        }
    }

    # Resetting non-empty default values for empty arrays and hashes
    my $create =  exists $self->_state->{CREATE} ? $self->_state->{CREATE} : undef;
    foreach my $option (@options) {
        my $type = $self->get_type($option);
        my $ref  = $self->get($option);
        if ($type eq 'hash') {
            if (defined $create && $option =~ /$create/) { # AppConfig bug: self-created hashes are polluted by key "1"
                delete $ref->{"1"} if (exists $ref->{"1"} && !defined $ref->{"1"});
            }
            $self->_default($option) if (!scalar keys %{ $self->get($option) });
        } elsif ($type eq 'array') {
            shift(@$ref) if (defined $create && $option =~ /$create/ && $ref->[0] eq "1"); # AppConfig bug: self-created arrays are polluted
            $self->_default($option) if (!scalar @{ $self->get($option) });
        }
    }

    # Put "-" back
    #@ARGV = map { s/_AppConfig_IGNORE_UNKNOWN_/-/g; $_ } @ARGV;

    # Expand @ARGV and STDIN if requested
    &_expand_array(defined $args && scalar @$args ? $args : \@ARGV) if ($self->{"EXPAND_ARGV"});
    &_expand_stdin() if ($self->{"EXPAND_STDIN"});

    # Show configuration
    if (grep { $_ eq 'debug' } @options) {
        $self->_dump_to_stderr() if ($self->debug =~ /^config/);
    }

    # Save configuration
    if ($self->configdump) {
        my @selected = grep { !/^(configfile|debug|help|configdump)$/ } $self->get_options;
        print $self->to_string(@selected);
        exit 0;
    }

    # Print help messages, if requested
    #  why did I use this test before, instead of just -t STDIN ????
    #	 (-t STDIN || (eof STDIN && !$self->{'EXPAND_STDIN'}))
    $self->{"EXPECT_INPUT"} = $self->_state->{EXPECT_INPUT} if (exists $self->_state->{EXPECT_INPUT});
    if ($self->help || ($self->{'EXPECT_INPUT'} && !(defined $args ? scalar(@$args) : scalar(@ARGV)) && -t STDIN)) {
        if (scalar @message) {
            $self->_pod(0, @message);
            return 1;
        } else {
            $self->_pod(1);
            exit 1;
        }
    }

    return 1;
}

#============================================================================
# Based on AppConfig::Getopt.pm by Andy Wardley <abw@wardley.org>
#============================================================================

=head2 getargs

 Title   : getargs
 Usage   : $config->getargs
 Function: parse arbitrarily ordered command line arguments
 Returns :
 Args    : array to process (default is @ARGV)

 After parsing, all options and parameters created using
 L<define> are removed from ARGV.

=cut

sub getargs {
    my $self  = shift; #
    my $args  = shift; # next parameter may be a reference to a list of args
    my $state = $self->_state;

    if (defined $args) {
        $self->{'_last_cmdline'} = [ $0, @$args ];
    } else {
        $self->{'_last_cmdline'} = [ $0, @ARGV ];
    }

    local $" = ', ';

    # we trap $SIG{__WARN__} errors and patch them into AppConfig::State
    local $SIG{__WARN__} = sub {
        my $msg = shift;
        # AppConfig::State doesn't expect CR terminated error messages
        # and it uses printf, so we protect any embedded '%' chars
        chomp($msg);
        #	$msg = "******* ERROR DETECTED WHILE PARSING COMMAND LINE OPTIONS *******\n\n$msg\n\n";
        #       $self->_pod(1, undef, $msg);
        die "$msg\n";
        exit 1;
    };

    # slurp all config items into @config
    my @config = ("no_auto_abbrev"); # Our default
    push(@config, "pass_through") if $self->{IGNORE_UNKNOWN};
    push(@config, "no_ignore_case") if $self->_state->{CASE};
    push(@config, shift) while defined $_[0] && ! ref($_[0]);

    # add debug status if appropriate (hmm...can't decide about this)
    #push(@config, 'debug') if $state->_debug();

    my @tmp = ();
    if (defined $args && $args != \@ARGV) {
        # copy any args explicitly specified into @ARGV
        @tmp = @ARGV;
        @ARGV = @$args;
    }

    # we enclose in an eval block because constructor may die()
    my ($getopt);
    eval {
        # configure Getopt::Long
        Getopt::Long::Configure(@config);

        # construct options list from AppConfig::State variables
        my @opts = $self->_getopt_state();

        # DEBUG
        if ($state->_debug()) {
            print STDERR "Calling GetOptions(@opts)\n";
            print STDERR "\@ARGV = (@ARGV)\n";
        };

        # call GetOptions() with specifications constructed from the state
        $getopt = GetOptions(@opts);
    };
    if ($@) {
        chomp($@);
        $state->_error("%s", $@);
        return 0;
    }

    # udpdate any args reference passed to include only that which is left
    # in @ARGV
    if (defined $args && $args != \@ARGV) {
        @$args = @ARGV;
        @ARGV  = @tmp;
    }

    return $getopt;
}

#------------------------------------------------------------------------
# Based on AppConfig::Getopt.pm by Andy Wardley <abw@wardley.org>
# _getopt_state()
#
# Constructs option specs in the Getopt::Long format for each variable
# definition.
#
# Returns a list of specification strings.
#------------------------------------------------------------------------

sub _getopt_state {
    my $self = shift;
    my ($var, $spec, $args, $argcount, @specs);

    my $linkage = sub { $_[1] = "$_[1]=$_[2]" if (scalar(@_) == 3); $self->_state->set(@_) };

    foreach $var (keys %{ $self->_state->{ VARIABLE } }) {
        $spec  = join('|', $var, $self->get_aliases($var));

        # an ARGS value is used, if specified
        unless (defined ($args = $self->_state->{ ARGS }->{ $var })) {
            # otherwise, construct a basic one from ARGCOUNT
            ARGCOUNT: {
                last ARGCOUNT unless defined ($argcount = $self->_state->_argcount($var));
                $args = "=s",  last ARGCOUNT if $argcount eq ARGCOUNT_ONE;
                $args = "=s@", last ARGCOUNT if $argcount eq ARGCOUNT_LIST;
                $args = "=s%", last ARGCOUNT if $argcount eq ARGCOUNT_HASH;
                $args = "!";
            }
        }
        $spec .= $args if defined $args;

        push(@specs, $spec, $linkage);
    }

    return @specs;
}

=head2 Internal methods

=head2 _dump_to_stderr

 Title   : _dump_to_stderr
 Usage   : $config = $config->_dump_to_stderr()
 Function: print a dump of the configuration
 Returns : nothing
 Args    :

 This method implicates the following options are added
 to your application interface:

 Option    : --debug
 Aliases   :
 Arguments :
 Default   : off

=cut

#------------------------------------------------------------------------
# _dump_to_stderr()
#
# Dumps the contents of the Config object and all stored variables.
#------------------------------------------------------------------------

sub _dump_to_stderr {
    my $self = shift;

    print STDERR "=======================================================================
    # AppConfig::AutoDoc instance INTERNAL STATE
    #
    # This class inherits most of its methods from Appconfig
    #
    # You may use this dump to inspect the parameters given to your program
    =======================================================================

    ";

    print STDERR "=" x 71, "\n";
    print STDERR
    "Status of AppConfig::State (version $VERSION) object:\n\t$self\n";


    print STDERR "- " x 36, "\nINTERNAL STATE:\n";
    foreach (qw( CREATE CASE PEDANTIC EHANDLER ERROR ), keys %EXTENDED_CLASS_TAGS) {
        printf STDERR "    %-14s => %s\n", $_,
        defined($self->_state->{ $_ }) ? $self->_state->{ $_ } : defined($self->{ $_ }) ? $self->{ $_ } : "<undef>";
    }

    print STDERR "- " x 36, "\nVARIABLES:\n";
    foreach my $var ($self->get_options) {
        $self->_dump_var_to_stderr($var);
    }

    print STDERR "- " x 36, "\n", "ALIASES:\n";
    my %aliases = $self->get_aliases;
    foreach my $var (sort keys %aliases) {
        printf STDERR "    %-12s => %s\n", $var, join(", ",@{ $aliases{$var} });
    }

    print STDERR "- " x 36, "\n", "ARGV: ",join(" ",@ARGV),"\n";
    print STDERR "=" x 72, "\n";

    exit 1;
}

=head2 _options_description

 Title   : _options_description
 Usage   : $config->_options_description
 Function: describes each option
 Returns : nothing
 Args    : none

 This method implicates the following options are added
 to your application interface:

 Option    : --usage
 Aliases   : --help --doc -h
 Arguments : none
 Default   : off

=cut

sub _options_description {
    my $self = shift;

    # Header for options
    my $text = "\n=head1 OPTIONS\n\nThis section describes program options, aliases and default values\n\n=head2 Program options\n\n";

    # Add options
    my $create = exists $self->_state->{'CREATE'} ? $self->_state->{'CREATE'} : undef;
    foreach my $option (grep { !$self->_private($_) } $self->get_options) {
        next if (defined $create && $option =~ /$create/); # Do not include auto-generated options
        $text .= "\n=head3 ".join(", ", map { length($_) == 1 ? "B<-$_>" : "B<--$_>" } ($option, sort { length($b) <=> length($a) } $self->get_aliases($option)))."\n";
        my $summary = $self->_extended_attribute($option,'SUMMARY') || " ";
        $text .= "\n$summary\n";
        $text .= "\n=over\n";
        $text .= "\n=item default: "._defaults2string($self,$option)."\n";
        $text .= "\n=back\n";
    }

    # Load predefined options
    $text .= "\n=head2 Options added by the configuration parser\n\nOptions below are generated by AppConfig::AutoDoc\n\n";
    foreach my $option (grep { $self->_private($_) } $self->get_options) {
        my $summary = $self->_extended_attribute($option,'SUMMARY') || " ";
        $text .= "\n=head3 ".join(", ", map { length($_) == 1 ? "B<-$_>" : "B<--$_>" } ($option, sort { length($b) <=> length($a) } $self->get_aliases($option)))."\n";
        $text .= "\n$summary\n";
        $text .= "\n=over\n";
        $text .= "\n=item default: "._defaults2string($self,$option)."\n";
        $text .= "\n=back\n";
    }

    return $text;
}


=head2 _options_summary

 Title   : _options_summary
 Usage   : $config->_options_summary
 Function: print application options summary
 Returns : nothing
 Args    : none

 This method implicates the following options are added
 to your application interface:

 Option    : --usage
 Aliases   : --help --doc -h
 Arguments : none
 Default   : off

=cut

sub _options_summary {
    my $self = shift;

    # AppConfig::State accessor
    my $create = exists $self->_state->{'CREATE'} ? $self->_state->{'CREATE'} : undef;

    # Load options data
    my @message = ([ "Long name","Aliases","Type",-1 ]);
    my @maxLength = map { length($_) } @{$message[0]};
    foreach my $option ($self->get_options) {
        # Do not include auto-generated options
        next if (defined $create && $option =~ /$create/);

        # Prepare and add this row
        my $aliases = join(", ",map { length($_) > 1 ? "--$_" : "-$_" } $self->get_aliases($option));
        my $type    = $self->get_type($option) || 'unknown';
        my $private = $self->_private($option);
        $option = length($option) > 1 ? "--$option" : "-$option";
        my $row = [ $option, $aliases, $type, $private ];
        foreach my $col (0..$#{$row}) { # Recalculate maximum column length
            $maxLength[$col] = length($row->[$col]) > $maxLength[$col] ? length($row->[$col]) : $maxLength[$col];
        }
        push(@message, $row);
    }

    # Right align each column and print row
    my $last = -10;
    my $ret .= "\n\n=head2 Program options summary\n\n";
    foreach my $row (sort { $a->[$#{$a}] <=> $b->[$#{$b}] || $a->[0] cmp $b->[0] } @message) {
        my @row = ();
        for (my $i=0; $i<=$#{$row}-1; $i++) {
            my $content = defined $row->[$i] ? $row->[$i] : '';
            push(@row, sprintf("%$maxLength[$i]s",$content));
        }
        $ret .=  "\n ".join(" : ", map { sprintf("%$maxLength[$_]s"," ") } 0..($#{$row}-1)) if ($last > -2 && $last != $row->[$#{$row}]);
        $ret .=  "\n ".join(" : ",@row);
        $last = $row->[$#{$row}];
    }
    $ret .= "\n";

    # Terminate the program
    return $ret;
}

=head2 _pod

 Title   : _pod
 Usage   : $config = $config->_pod(1)
 Function: show detailed program documentation
 Returns : nothing
 Args    : (boolean) default 1, (de)activate pager
           (string) optional, format this POD text
                    instead of the caller program ($0)
           (array of strings) optional, append

 This method parses any POD documentation added to a
 calling program and adds descriptions of the program
 options

 This method implicates the following options are added
 to your application interface:

 Alternatively, you may either format or append any number
 of POD to the formatted documention printed

 Option    : --doc
 Aliases   :
 Arguments : none
 Default   : off

=cut

sub _pod {
    my ($self, $pager, $text, @append) = @_;

    # open this file
    unless (defined $text) {
        $text = "";
        if ( -e $0 ) {
            open(F,"<$0") || die "Could not generate documentation for $0:\nPerhaps you don't have authorization to read this program's code (read flag off)?";
            $text = join("", <F>);
            close(F);
        }

        # Include summary and extra mesage at the bottom
        $text .= $self->_options_description;
        $text .= $self->_options_summary;
    }

    # Append something else
    foreach my $append (@append) {
        $text .= $append if (defined $append);
    }

    # Redirect STDOUT
    if ($pager) {
        if (-t STDOUT) {
            if (exists $ENV{"PAGER"}) {
                $ENV{'PAGER'} = "less -r" if ($ENV{'PAGER'} =~ m|/?less$|);
                if (open(PAGER,"| $ENV{PAGER}")) {
                    *STDOUT = *PAGER;
                }
            } elsif (open(PAGER,"| less -r")) {
                *STDOUT = *PAGER;
            } elsif (open(PAGER,"| more")) {
                *STDOUT = *PAGER;
            }
        } elsif (-t STDERR) {
            my $OriginalCmdLine = join(" ",@{ $self->{'_last_cmdline'} });
            $self->_default_error_handler("Missing or inapropriate arguments triggered an attempt to dump the program's documentation to a file or pipe:\n\n$OriginalCmdLine\n\nTry '$0 --help' to learn how to use it\nAbort!!!");
        }
    }

    # Create IO::File handler for Pod::Text::Termcap
    require IO::String;
    my $io = IO::String->new($text);

    # Generate actual documentation
    require Pod::Text::Termcap;
    my $parser = Pod::Text::Termcap->new(sentence => 0, width => 78);
    #my $parser = Pod::Man->new();
    $parser->parse_from_filehandle($io);

    # Terminate
    close(PAGER) if ($pager);
    return 1;
}

sub _defaults2string {
    my ($self,$option) = @_;

    # Laod defaults
    my $defaults = $self->get_defaults($option);

    # Convert default lists to strings
    if (ref($defaults) eq "ARRAY") {
        $defaults = scalar(@{$defaults}) ? join(", ",@{$defaults}) : '<empty>';
    } elsif (ref($defaults) eq "HASH") {
        $defaults = scalar(keys %{$defaults}) ? join(", ",map { $_."=".$defaults->{$_} } keys %{$defaults}) : '<empty>';
    } elsif (!$self->_argcount($option) && !$defaults) {
        $defaults = "off";
    }

    # Check defaults for non-printable characters
    if (!defined $defaults) {
        $defaults = "<undef>";
    } elsif ($defaults eq "\t") {
        $defaults = "<tab>";
    }

    return $defaults;
}

=head2 _state

 Title   : _state
 Usage   : my $state = $config->_state
 Function: get the internal AppConfig::State object
 Returns : AppConfig::State
 Args    :

=cut

sub _state {
    return shift->{STATE};
}

=head2 _no_expect_input

 Title   : _no_expect_input
 Usage   :
 Function: turn off input check
 Returns :
 Args    :

=cut

sub _no_expect_input {
    my ($state,$name,$value) = @_;
    $state->{'EXPECT_INPUT'} = 0;
}

=head2 _expand_file

 Title   : _expand_file
 Usage   :
 Function: Internal callback function for EXPAND_FILES tag
           Replaces filenames with their contents
 Returns : true
 Args    :

=cut

sub _expand_file {
    my ($state,$name,$value) = @_;
    my $ref = $state->get($name);
    my $ret = undef;
    my $type = ref($ref);
    if ($type eq 'ARRAY') {
        $ret = _expand_array($ref,$#{$ref});
    } elsif ($type eq 'HASH') {
        $ret = _expand_hash($ref,"file");
    } elsif ($state->_argcount($name) == ARGCOUNT_ONE) {
        return 0 unless (defined $value && (-e $value || $value eq '-'));
        open(SCALAR,"<$value") || die "Could not expand file $value";
        my $string = join("",<SCALAR>);
        close(SCALAR);
        # Once more I'm forced to bypass AppConfig's interface due to
        # lack of functionality. In this case, we have no way to set
        # the value of scalars without running
        $state->{ VARIABLE }->{ $name } = $string;
        #$state->set($name,$string);
    }
    return $ret;
}

=head2 _expand_array

 Title   : _expand_array
 Usage   : $config->_expand_array
 Function: replaces filenames in @ARGV by their contents
 Returns : true
 Args    :

=cut

sub _expand_array {
    my ($ref,$elm) = @_;
    my ($i,$max) = defined($elm) ? ($elm,$elm) : (0,$#{$ref});
    while ($i <= $max) {
        my $file = $ref->[$i];
        if (-f $file || $file eq '-') {
            return 0 if ($file eq '-' && -t STDIN);
            open(FILE,"<$file") || die "Could not open file $file";
            my @add = map { chomp; s/^\s+//; s/\s+$//; $_ } <FILE>;
            close(FILE);
            splice(@$ref,$i,1,@add);        # Replace element with file content
            $i += scalar(@add);             # Next or last (if $#add is still 0)
            $max = $#{$ref} unless (defined $elm);
        } else {
            $i++;
        }
    }
    return $i;
}

=head2 _expand_hash

 Title   : _expand_hash
 Usage   : $config->_expand_hash
 Function: loads a two column file into a hash
 Returns : true
 Args    :

=cut

sub _expand_hash {
    my ($ref,@keys) = @_;
    @keys = keys %$ref unless (scalar @keys);
    my $ret = 0;
    foreach my $key (@keys) {
        next unless (exists $ref->{$key});
        my $file = $ref->{$key};
        die "Unable to load hash from file $file" if (! -e $file || ($file eq '-' && -t STDIN));
        delete $ref->{$key};
        open(FILE,"<$file") || die "Could not open file $file";
        while (<FILE>) {
            chomp; s/^\s+//; s/\s+$//; # Remove newline, leading and trailing spaces
            my ($key,$value) = split(/$EXPAND_FILES_OPTIONS{COLUMN_DELIMITER}/,$_,2);
            $ref->{$key} = $value;
            $ret++;
        }
        close(FILE);
    }
    return $ret;
}

=head2 _expand_stdin

 Title   : _expand_stdin
 Usage   : $config->_expand_stdin
 Function: replaces filenames in @ARGV by their contents
 Returns : true
 Args    :

=cut

sub _expand_stdin {
    return 0 if (-t STDIN);
    my @tmp = @ARGV;
    @ARGV = ();
    @ARGV = map { chomp; s/^\s+//; s/\s+$//; $_ } <>;
    push(@ARGV, @tmp);
}

=head2 _extended_attribute

 Title   : _extended_attribute
 Usage   : $slef->_extended_attribute("help",'SUMMARY',"print help message");
 Function: get/set extended option attributes
 Returns : string or undef
 Args    : option name, tag name and (optional) new value

=cut

sub _extended_attribute {
    my ($self, $option, $tag, $value) = @_;
    if (defined $value) {
        if ($tag eq 'EXPAND_FILES') {
            $self->{'_attributes'}{$tag}{$option} = $value;
        } else {
            $self->{'_attributes'}{$tag}{$option} = $value;
        }
    } elsif (!defined $self->{'_attributes'}{$tag}{$option}) {
        $self->{'_attributes'}{$tag}{$option} = $EXTENDED_OPTION_ATTRIBUTES{$option};
    }
    return $self->{'_attributes'}{$tag}{$option};
}

=head2 _private

 Title   : _private
 Usage   : $slef->_private("help");
 Function: get/set option as private to this class
 Returns : boolean
 Args    : none

=cut

sub _private {
    my ($self, $option) = @_;
    return exists($DEFAULT_OPTIONS{$option});
}

=head2 _config2hash

 Title   : _config2hash
 Usage   : %hash = $self->_config2hash(1);
 Function: copies the values of all configuration
           parameters in a hash, and return it
 Returns : hash
 Args    : (boolean) when true, tells the method to
           remove default values from hashes and arrays,
           thus allowing fresh additional rounds of
           parsing

=cut

sub _config2hash {
    my ($self, $clean, @keys) = @_;
    $clean = 0 unless (defined $clean);

    # Store all non-empty value from a configuration object
    my %hash = ('_defaults' => []);
    my %config = $self->varlist('.');
    @keys  = keys %config unless (scalar @keys); # process all vars unless set by caller
    foreach my $option (@keys) {
        #next if ($self->_private($option)); # Do not include private options

        # Load option information
        my $value    = $config{$option};
        my $defaults = $self->get_defaults($option);
        my $type     = $self->get_type($option);

        # If $defaults and/or value is undefined
        if (!defined $value) {
            if (!defined $defaults) {
                push(@{ $hash{"_defaults"} }, $option);
            }
            next;
        }

        # Delete options set to default values and copy the rest
        if ($type eq 'boolean' || $type eq 'scalar') {
            if (defined($defaults) && $value eq $defaults) {
                push(@{ $hash{"_defaults"} }, $option);
            } else {
                $hash{$option} = $value;
            }
        } elsif ($type eq 'array') {
            $defaults = [] if (!defined $defaults);
            $defaults = join("\cA",sort map { $_ || "" } @{$defaults});
            my $test  = join("\cA",sort map { $_ || "" } @{$value});
            if ($test eq $defaults) {
                push(@{ $hash{"_defaults"} }, $option);
                @{ $self->get($option) } = () if ($clean);
            } elsif (scalar(@$value)) {
                $hash{$option} = [ @{$value} ];
            }
        }
        elsif ($type eq 'hash') {
            $defaults = {} if (!defined $defaults);
            $defaults = join("\cA",map { $_."\cB".$defaults->{$_} } sort keys %{$defaults});
            my $test  = join("\cA",map { $_."\cB".$value->{$_} } sort keys %{$value});
            if ($test eq $defaults) {
                push(@{ $hash{"_defaults"} }, $option);
                %{ $self->get($option) } = () if ($clean); # Cleaning deafult values from hash
            } elsif (scalar(keys %$value)) {
                $hash{$option} = { %{$value} };
            }
        } else {
            die "Unknown type $value associated with option $option";
        }
    }

    return %hash;
}

=head2 _reset_from_hash

 Title   : _reset_from_hash
 Usage   : %hash = $slef->_reset_from_hash();
 Function: reset option values based on a hash created
           by _config2hash
 Returns :
 Args    : (boolean) when true, resets default values for
           hashes that were previouly
           set to default values  reported

 Note: used with _config2hash to reset user
       defined options to a previous state

=cut

sub _reset_from_hash {
    my ($self,%hash) = @_;

    foreach my $option (grep { $_ ne '_defaults' } keys %hash) {
        $self->_default($option);
        RESET: {
            $_ = $self->get_type($option);
            /scalar|boolean/ && do {
                $self->set($option, $hash{$option});
                last RESET;
            };
            /array/ && do {
                @{ $self->get($option) } = @{$hash{$option}};
                last RESET;
            };
            /hash/ && do {
                %{ $self->get($option) } = %{$hash{$option}};
                last RESET;
            };
            die "Unknown type $_ associated with option $option";
        }
    }
}

=head2 _dump_var_to_stderr

 Title   : _dump_var_to_stderr
 Usage   : $self->_dump_var_to_stderr();
 Function: displays the content of the specified variable
 Returns :
 Args    : a variable name or aliases

 Note: this is a modified version of _dump_var from AppConfig::State

=cut

sub _dump_var_to_stderr {
    my $self   = shift;
    my $var    = shift;

    return unless defined $var;

    # $var may be an alias, so we resolve the real variable name
    my $real = $self->_state->_varname($var);
    if ($var eq $real) {
        print STDERR "$var\n";
    }
    else {
        print STDERR "$real  ('$var' is an alias)\n";
        $var = $real;
    }

    # for some bizarre reason, the variable VALUE is stored in VARIABLE
    # (it made sense at some point in time)
    printf STDERR "    VALUE        => %s\n", defined($self->_state->{ VARIABLE }->{ $var })
    ? ref($self->_state->{ VARIABLE }{ $var }) eq "HASH"
    ? "( ".join(", ",map { $_." => ".$self->_state->{ VARIABLE }{ $var }{ $_ } } sort keys %{ $self->_state->{ VARIABLE }{ $var } })." )"
    : ref($self->_state->{ VARIABLE }{ $var }) eq "ARRAY"
    ? "( ".join(", ", @{ $self->_state->{ VARIABLE }{ $var } })." )"
    : $self->_state->{ VARIABLE }->{ $var }
    : "<undef>";

    # the rest of the values can be read straight out of their hashes
    foreach my $param (qw( DEFAULT ARGCOUNT VALIDATE ACTION EXPAND ), sort keys %EXTENDED_OPTION_ATTRIBUTES) {
        my $ref = exists $EXTENDED_OPTION_ATTRIBUTES{$param} ? $self->{'_attributes'} : $self->_state;
        printf STDERR "    %-12s => %s\n", $param, defined($ref->{ $param }->{ $var })
        ? ref($ref->{ $param }{ $var }) eq "HASH"
        ? "( ".join(", ",map { $_." => ".$ref->{ $param }{ $var }{ $_ } } sort keys %{ $ref->{ $param }{ $var } })." )"
        : ref($ref->{ $param }{ $var }) eq "ARRAY"
        ? "( ".join(", ", @{ $ref->{ $param }{ $var } })." )"
        : $ref->{ $param }{ $var }
        : "<undef>";
    }

    # summarise all known aliases for this variable
    print STDERR "    ALIASES      => ",join(", ",$self->get_aliases($var)),"\n";
}

=head2 _default_error_handler

 Title   : _default_error_handler
 Usage   : new({ ERROR => \&_default_error_handler })
 Function:
 Returns : nothing
 Args    : none

=cut

sub _default_error_handler {
    my $format  = shift;
    my $message = _stack_trace_dump(@_);
    die $message;
}

=head2 Methods borrowed from bioperl (Bio::Root::RootI)

=head2 _stack_trace_dump

 Title   : _stack_trace_dump
 Usage   :
 Function:
 Example :
 Returns :
 Args    :


=cut

sub _stack_trace_dump {
    my @message = @_;
    my @stack = _stack_trace();

    shift @stack;
    shift @stack;
    shift @stack;

    my $out;
    my ($module,$function,$file,$position);
    foreach my $stack ( @stack) {
        ($module,$file,$position,$function) = @{$stack};
        $out .= "** STACK $function $file:$position\n";
    }
    $out .= "\n";
    $out .= join("\n",@message)."\n";

    return $out;
}

=head2 _stack_trace

 Title   : _stack_trace
 Usage   : @stack_array_ref = _stack_trace()
 Function: gives an array to a reference of arrays with stack trace info
           each coming from the caller(_stack_number) call
 Returns : array containing a reference of arrays
 Args    : none


=cut

sub _stack_trace {
    my $i = 0;
    my @out = ();
    my $prev = [];
    while( my @call = caller($i++)) {
        # major annoyance that caller puts caller context as
        # function name. Hence some monkeying around...
        $prev->[3] = $call[3];
        push(@out,$prev);
        $prev = \@call;
    }
    $prev->[3] = 'toplevel';
    push(@out,$prev);
    return @out;
}

# Return true to require and use builtin perl functions
1;
