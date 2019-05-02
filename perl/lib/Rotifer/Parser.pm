=head1 NAME

Rotifer::Parser - parser factory for Rotifer 

=head2 SYNOPSIS

 my $parser = Rotifer::Parser->create("arch");
 my @data_stream = $parser->parse(@ARGV);
 foreach my $data (@data_stream) {
   ... do something ...
 }

=head2 DESCRIPTION

Rotifer::Parser is a factory that to load modules that implement
input file parsers for many Rotifer components.

This factory uses MooseX::AbstractFactory for loading modules and
enforces a very basic interface for parsers (see below). 

=head2 PARSER METHODS

This factory tests the interface of parsers for the implementation
of the following methods: 

=head3 parse

This method should accept as input an array of filenames or file
handles to process. Optionally, a reference to a hash of options
may also be accepted as a first argument but the interface of the
parser has to document this feature.

It is expected to return an array whose contents are defined by
the parser. If no valid input is found, the parser should return
an empty array and throw an exception.

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Rotifer::Parser;

#use MooseX::AbstractFactory;
use Module::Load; # 'load' method
use Carp::Clan qr/^Rotifer::Parser/;
use Rotifer::Utils qw(describe_subclasses);
use strict;
use warnings;

=head2 METHODS

=head2 create

 Title   : create
 Usage   : Rotifer::Parser->create("arch")
 Function: load and maybe return a parser module
 Returns : whatever the module's new() method returns,
           if it exists. Undef ortherwise. 
 Args    : arguments to the parser's "new" method

=cut

sub create {
    my ($self,$parser,@args) = @_;
    confess "Unable to guess parser!" if (!defined $parser);
    if ($parser eq "help") {
	print describe_subclasses(__PACKAGE__);
	exit 1;
    }
    $parser = "${self}::$parser";
    load "$parser" || confess "Unable to load $parser!";
    return $parser->new(@args);
}

# Set my subclass namespace
#implementation_class_via sub { 
#    my $parser = shift;
#    if ($parser eq "help") {
#	print describe_subclasses(__PACKAGE__);
#	exit 1;
#    }
#    return __PACKAGE__."::$parser";
#};

#__PACKAGE__->meta->make_immutable;
1;
