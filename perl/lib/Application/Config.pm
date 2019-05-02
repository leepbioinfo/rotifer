# POD documentation - main docs before the code

=head1 NAME

Application::Config - parse command line options and/or configuration files

=head1 SYNOPSIS

  # Simple example for single value (scalar) option
  use Application::Config qw(:argcount GetConfig);
  my $config = GetConfig('opt1' => {
      ALIAS    => "o",
      ARGCOUNT => ARGCOUNT_ONE,
      DEFAULT  => "<undef>",
      SUMMARY  => "simple scalar option",
  });

 # Accessing values in the program
 my $parameter = $config->opt1;

 # Running the program
 my_prog --opt1 20

 See AppConfig::AutoDoc and AppConfig for additional examples

=head1 DESCRIPTION

This module standadizes parsing of program options from the command
line and configuration file(s). It provides automatic formatting of
short documentation messages for each option and parsing of the
program's embedded POD documentation, if availbale.

See the documentation below for information on how to create and 
define your application's options and read the documentation of the
B<AppConfig> and B<AppConfig::AutoDoc> libraries for more details on
the parameters accepted by GetConfig.

=head2 Documentation

Application::Config was designed as wrapper to build AppConfig::AutoDoc
objects that do the actual command line and file parsing. Therefore, 
its documentation is not complete, as most of the information is in the
libraries it depends on. In order to become a power user, you might read
or search for the features you need in the documentation for
Application::Config's depedencies. The dependency tree is:

 Application::Config
  AppConfig::AutoDoc
   AppConfig
    AppConfig::State
    AppConfig::File

=head2 Precedence

Options given by users to an application based in this module have
the following precedence:

1) command line 

2) configuration files

3) program hard-coded defaults

If a default configuration is set in the program's code, options 
loaded from alternative configuration files choosen using option
--configfile (see below) are given precedence.

=head2 Options defined by this module

The following set of predefined options is provided by this module:

 --usage or --help or -h
 --doc
 --debug
 --configfile

Options --usage and --doc trigger execution of internal routines
that automatically generate and print documentation added by the
user to its program in POD format and derived by AppConfig::AutoDoc
from the program option's definitions.

Option --configfile directs the parser to process one or more 
configuration files. Parsing of these files is automatic and
more than one configuration file may be chosen by repeating 
this option in the command line (--configfile file1 --configfile
file2).

Option --debug does not trigger execution of any predefined method
unless its argument is 'config' (--debug config). When 'config' is
given as argument to --debug, all information about the configured
parameters is formatted and printed to STDERR.

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Application::Config;

use AppConfig::AutoDoc qw(:argcount :expand);
use Exporter;
use base qw(Exporter);
use strict;
use warnings;

our $VERSION = "0.1";

# Expanding AppConfig's exported symbols list
our @EXPORT_OK   = (@AppConfig::EXPORT_OK);
our %EXPORT_TAGS = (%AppConfig::EXPORT_TAGS);
push(@EXPORT_OK, qw(GetConfig GetTableToolsConfig expand_numeric_interval quote_string));

=head2 GetConfig

 Title   : GetConfig
 Usage   : $config = GetConfig()
 Function: parse input and build a configuration object
 Returns : AppConfig::AutoDoc object
 Args    : hash reference (optional) and/or
           hash of hashes defining options (optional)

 Note: This method prepares a configuration object of the
       class AppConfig::AutoDoc and uses it to parse the 
       command line and configuration files set by
       --configfile. It accepts the same set of arguments
       as AppConfig::AutoDoc::new. 

       See the documentation for the methods L<new> and 
       L<define> in AppConfig for a detailed description of
       how to define options.

=cut

# This code is a modified version of Bio::SeqIO::new
sub GetConfig {
    my (@args) = @_;
    my $self = shift(@args) if (ref($args[0]) =~ /^\S+\=\S+\([\dx]\)/);
    my $config = AppConfig::AutoDoc->new(@args);
    $config->parse_cmdline_files();
    return $config;
}

=head2 GetTableToolsConfig

 Title   : GetTableToolsConfig
 Usage   : $config = GetTableToolsConfig()
 Function: parse input and build a configuration object
 Returns : AppConfig::AutoDoc object
 Args    : hash reference (optional) and/or
           hash of hashes defining options (optional)

 Note: This method prepares a configuration object of the
       class AppConfig::AutoDoc and uses it to parse the 
       command line and configuration files set by
       --configfile. It accepts the same set of arguments
       as AppConfig::AutoDoc::new. 

       See the documentation for the methods L<new> and 
       L<define> in AppConfig for a detailed description of
       how to define options.

       This function extends the base configuration object
       created by GeConfig (above) with the common set of
       options used by Rotifer's table tools like tjoin, 
       tgroup, etc. The following options are added:

       Option    : --empty
       Aliases   : -e
       Arguments : String to add to empty cells
       Default   : off

       Option    : --input_delimiter
       Aliases   : -s
       Arguments : string or regular expression that delimits
                   (separates) input fields (i.e. columns)
       Default   : <tab> (\t)

       Option    : --output_delimiter
       Aliases   : -r
       Arguments : string or regular expression that delimits
                   (separates) input fields (i.e. columns)
       Default   : <tab> (\t)

       Option    : --pad
       Aliases   : -p
       Arguments : Add spaces so that output columns have equal lengths
       Default   : off

       Option    : --parse_header
       Aliases   : -y
       Arguments : Parse first row in each file as table header
       Default   : off

       Option    : --rename_duplicates
       Aliases   : -rd
       Arguments : Add <_number> to repeated column names
       Default   : off

       Option    : --sort
       Aliases   : -o
       Arguments : List of columns to sort the output table
       Default   : 0 (first column)

       Option    : --unpad
       Aliases   : -u
       Arguments : Remove trailing and leading spaces from input values
       Default   : on

=cut

# This code is a modified version of Bio::SeqIO::new
sub GetTableToolsConfig {
    my (@args) = @_;
    my $self = shift(@args) if (ref($args[0]) =~ /^\S+\=\S+\([\dx]\)/);
    my $opts = shift(@args) if (ref($args[0]) eq "HASH");

    @args = ( # Table options

	      'empty' => {
		  DEFAULT  => '',
		  ARGCOUNT => ARGCOUNT_ONE,
		  ALIAS    => 'e',
		  SUMMARY  => "string to fill empty cells",
	      },

	      'input_delimiter' => {
		  DEFAULT  => "\t",
		  ARGCOUNT => ARGCOUNT_ONE,
		  ALIAS    => 's',
		  SUMMARY  => "column delimiter of input files",
	      },

	      'output_delimiter' => {
		  DEFAULT  => "\t",
		  ACTION   => \&quote_string,
		  ARGCOUNT => ARGCOUNT_ONE,
		  ALIAS    => 'r',
		  SUMMARY  => "column delimiter for output table",
	      },

	      'pad' => {
		  DEFAULT  => 0,
		  ARGCOUNT => ARGCOUNT_NONE,
		  ALIAS    => 'p',
		  SUMMARY  => "Add trailing spaces all OUTPUT rows in all columns so that each column's length is constant.",
	      },

	      'parse_header' => {
		  DEFAULT  => 0,
		  ARGCOUNT => ARGCOUNT_NONE,
		  ALIAS    => 'y|header',
		  SUMMARY  => "Boolean flag that indicates first row contain column names. See also documentation of option --rename_duplicates.",
	      },

	      'rename_duplicates' => {
		  DEFAULT  => 0,
		  ACTION   => sub { my ($s,$n,$v) = @_; $s->set("header",1) },
		  ARGCOUNT => ARGCOUNT_NONE,
		  ALIAS    => 'rd',
		  SUMMARY  => "Rename duplicated columns (activates --header).",
	      },

	      'sort' => {
		  ALIAS    => 'o',
		  ACTION   => \&expand_numeric_interval,
		  ARGCOUNT => ARGCOUNT_LIST,
		  DEFAULT  => [ 0 ],
		  SUMMARY  => "list of columns used to order the output table",
	      },

	      'unpad' => {
		  DEFAULT  => 1,
		  ARGCOUNT => ARGCOUNT_NONE,
		  ALIAS    => 'u',
		  SUMMARY  => "Remove leading/trailing invisible characters from all INPUT columns.",
	      },

	      @args
	);
	unshift(@args,$opts) if (defined $opts);

    # Same as GetConfig
    my $config = AppConfig::AutoDoc->new(@args);
    $config->parse_cmdline_files();

    return $config;
}

=head2 expand_numeric_interval

 Title   : expand_numeric_interval
 Usage   : ACTION => \&Application::Config::expand_numeric_interval
 Function: convert numeric interval to array
 Returns : nothing
 Args    : none

=cut

sub expand_numeric_interval {
    my ($s,$n,$v) = @_;
    if ($v =~ /\d+\.\.\d+/) {
	$s = $s->get($n);
	pop(@$s);
	push(@$s, eval "$v");
    }
}

=head2 quote_strng

 Title   : quote_strng
 Usage   : ACTION => \&Application::Config::quote_strng
 Function: quote string
 Returns : nothing
 Args    : none

=cut

sub quote_string {
    my ($s,$n,$v) = @_;
    my $new = eval qq{ return "$v" };
    $s->set($n,$new) if ($v ne $new);
}

1;
# Moose
#no Moose;
#__PACKAGE__->meta->make_immutable;
