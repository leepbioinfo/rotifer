#!/usr/bin/env perl

# LIBRARIES and PRAGMAS: start
use FindBin;
use lib "$FindBin::RealBin/../perl/lib";
use strict;
use warnings;
# LIBRARIES and PRAGMAS: end

# Parse configuration/options
my $CONFIG = &parse_configuration();

# Your code could replace the lines below!!!!!!!
# See the examples below on how to access the command line
# options from inside the program

# MAIN PROGRAM: start

# Dump config to a file
if (defined $CONFIG->dump) {
    $CONFIG->to_file("test.conf",grep { !/^(dump|help|debug)$/ && defined $CONFIG->$_ } $CONFIG->get_options);
    exit 0;
}

# Acessing the value of booleans (ARGCOUNT_NONE) and scalars (ARGCOUNT_ONE)
print "Option0 is ",$CONFIG->option0,"\nOption1 is ",$CONFIG->option1,"\n";

# Accessing values in arrays (ARGCOUNT_LIST)
print "Option2 is ( ",join(", ",@{ $CONFIG->option2 })," )\n";

# Accessing values in arrays (ARGCOUNT_LIST)
print "Option3 is ( ",join(", ",map { $_." => ".$CONFIG->option3->{$_} } sort keys %{ $CONFIG->option3 })," )\n";

# Inspecting values in @ARGV after parsing
print 'After parsing, @ARGV contains ( ',join(", ",@ARGV)," )\n";

exit 0;
# MAIN PROGRAM: end

###############
# Subroutines

# Creating the configuration object
#
# The Application::Config library we are using here is just
# a wrapper that facilitates the use of AppConfig::AutoDoc.
#
# Because it depends on other libraries,Application::Config
# documentation is far from complete but contains the pointers
# to detailed documentation on its dependencies. The user should
# start reading "perldoc Application::Config" and them read the
# docs on each of its dependencies in the oerder below:
#
# Application::Config
#  AppConfig::AutoDoc
#   AppConfig
#    AppConfig::State
#
# The argument for new is an anonymous hash. This hash is
# used to control the bahaviour of AppConfig::AutoDoc. In
# This case, setting EXPAND_ARGV to 1 enables tranfers of
# all words in input files to @ARGV. EXPAND_STDIN does the
# same for text comming through pipelines (standard input)
#
# See perldoc AppConfig::AutoDoc for details
sub parse_configuration {
    use Application::Config qw(:argcount GetConfig);
    # Here I'm passing GetConfig a hash that serves both as an example and
    # prototype for the kind of argument expected by the method
    # "define" (see below) from AppConfig::AutoDoc.
    my $appconfig = GetConfig(# First, lets configure AppConfig::Autodoc by
			      # giving a hash reference with parameters
			      { 
				  CASE           => 1, 
				  EXPAND_ARGV    => 1, 
				  EXPAND_STDIN   => 1, 
				  IGNORE_UNKNOWN => 1, 
				  ERASE_DEFAULTS => 1, 
				  CREATE         => '^alias_',
				  GLOBAL         => { ARGCOUNT => ARGCOUNT_HASH },
			      },

			      # Now, some examples of options you can create
			      "option0" => {
				  ALIAS    => "a",
				  ARGCOUNT => ARGCOUNT_NONE,
				  DEFAULT  => 0,
				  SUMMARY  => "a flag that is either true (1) or false (0)",
			      },
			      "option1" => {
				  ACTION   => 'EXPAND_FILES',
				  ALIAS    => "b",
				  ARGCOUNT => ARGCOUNT_ONE,
				  DEFAULT  => undef,
				  SUMMARY  => "a string set by -b value",
			      },
			      "option2" => {
				  ACTION   => 'EXPAND_FILES',
				  ALIAS    => "c",
				  ARGCOUNT => ARGCOUNT_LIST,
				  DEFAULT  => [],
				  SUMMARY  => "an array of values (-c value1 -c value2...)",
			      },
			      'option3' => {
				  ACTION   => 'EXPAND_FILES',
				  ALIAS    => "d",
				  ARGCOUNT => ARGCOUNT_HASH,
				  DEFAULT  => { "autodoc" => "test" },
				  SUMMARY  => "a hash, i.e. pairs of value like -d key1=value1 -d key2=value2...",
			      },
			      'section_array' => {
				  ACTION   => 'EXPAND_FILES',
				  ALIAS    => "sa",
				  ARGCOUNT => ARGCOUNT_LIST,
				  DEFAULT  => [],
				  SUMMARY  => "an array in section...",
			      },
			      'section_hash' => {
				  ACTION   => 'EXPAND_FILES',
				  ALIAS    => "sh",
				  ARGCOUNT => ARGCOUNT_HASH,
				  DEFAULT  => {},
				  SUMMARY  => "a hash in section...",
			      },
			      'dump' => {
				  ARGCOUNT => ARGCOUNT_ONE,
				  DEFAULT  => undef,
			      },
			      "configfile" => {
				  DEFAULT => [ "example.conf" ] 
			      },
	);

    return $appconfig;
}


# POD: start
# Documentation (use option -doc or --doc to see it!)
#
# AppConfig::AutoDoc can automatically extract Plain Old Documentation
# from within the caller program's, add descriptions of the options
# created by L<define> and do some pretty formatting for output.
# 
# Note that POD may be added anywhere in your program, thus allowing 
# the documentation for a program to be left side by side with the
# function's definition.

=head1 NAME

 test_appconfig_autodoc.pl - 
        example and prototype for AppConfig::AutoDoc

=head1 SYNOPSIS

 Setting the boolean flag on (by default is off)

 test_appconfig_autodoc.pl -option0

 Making sure the boolean flag is off

 test_appconfig_autodoc.pl -nooption0

 Setting parameter "scalar" to a string

 test_appconfig_autodoc.pl -option1 "a string"

 Setting array "option2" to a list of strings

 test_appconfig_autodoc.pl -option2 string1 -option2 string2

 Setting hash "option3" to pairs of values

 test_appconfig_autodoc.pl -option3 key1=value1 -option3 key2=value2

=head1 DESCRIPTION

This program provides an example of how to use AppConfig::AutoDoc

Please consider reading the code of this program and using it as a 
template to build your own AppConfig::AutoDoc based programs. Also,
please read the module's documentation (perldoc AppConfig::AutoDoc)
and, since most functionality in AppConfig::AutoDoc is derived from
AppConfig, see also AppConfig documentation (perldoc AppConfig)

=head1 AUTHOR

 Robson Francisco de Souza

=cut
# POD: end
