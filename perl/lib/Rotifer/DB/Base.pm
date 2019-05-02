# POD documentation - main docs before the code

=head1 NAME

Rotifer::DB::Base - base methods to discover how to
                    access your databases

=head1 SYNOPSIS

  # Retrieving FASTA sequences
  use Rotifer::DB::Base;
  my @data = Rotifer::DB::Base->datadir();

=head1 DESCRIPTION

Rotifer::DB::Base implements methods to find local
databases paths.

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::DB::Base;

use strict;
use warnings;
use autodie qw(:all);
use Carp::Clan qr/^Rotifer::DB/;
use File::Spec;
use Rotifer;
use Sub::Exporter -setup => {
    exports => [qw/data_path/],
};
use vars qw(@DATADIRS);

=head2 datadir

 Title   : datadir
 Usage   : Rotifer::DB::Base->datadir()
 Function: list of database directories
 Returns : array
 Args    : new list

=cut

sub datadir {
    my $self = shift if (defined $_[0] && ($_[0] eq __PACKAGE__ || ref($_[0]) eq __PACKAGE__));
    if (scalar @_) {
	@DATADIRS = @_;
    } else {
#	@DATADIRS = grep { defined && -d } (Rotifer->userdir("data"), $ENV{ROTIFER_DATA}, Rotifer->sysdir("share"))
	@DATADIRS = grep { defined && -d } ($ENV{ROTIFER_DATA}, Rotifer->sysdir("share"))
	    unless (scalar @DATADIRS);
    }
    return @DATADIRS;
}

=head2 data_path

 Title   : data_path
 Usage   : Rotifer::DB::Base->data_path("fadb","euk_fa")
 Function: find the first valid path to your local databases
 Returns : array, string or undef (context dependent)

  - scalar context : first directory where query was found
                     or undef if not found

  - array context  : all paths that match your query

 Args    : array of path components (strings)
 Note    : Search order is given by Rotifer::DB::Base->datadir()

=cut

sub data_path {
    my $self = shift if (defined $_[0] && ($_[0] eq __PACKAGE__ || UNIVERSAL::isa($_[0],__PACKAGE__)));

    my @found = ();
    my $path = File::Spec->catfile(@_);
    foreach my $dir (&datadir()) {
	my $fullpath = File::Spec->catfile($dir,$path);
	if ( -e $fullpath ) {
	    if (wantarray) {
		push(@found, $fullpath);
	    } else {
		return $fullpath;
	    }
	}
    }

    return wantarray ? @found : undef;
}

# Make perl compiler happy
1;
