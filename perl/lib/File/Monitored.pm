# POD documentation - main docs before the code

=head1 NAME

File::Monitored - create and remove monitored files

=head1 SYNOPSIS

  # Creating files that are removed if the program is terminated
  use File::Monitored qw(topen);
  my $fh = topen("data.txt");
  print $fh "Write nice things to data.txt";
  tclose($fh);

=head1 DESCRIPTION

Vanilla perl functions for creating and erasing file (close and
open) tend to leave files on disk if a program is interrupted
before finishing execution. Although this is ok in most cases ,
certain tasks may require creation of temporary or intermediate
files that are expected to be removed when the program finishes
and will just pollute the working directory if the program is
interrupted.

File::Monitored provides some basic subroutines to create files
that will be removed if the program is interrupted.

=head1 METHODS

This module does not exports anything by default, therefore use

use File::Monitored qw(something)

to have a method or set of methods called "something" made 
available. The following tags can be used to import several
at once into your program:

=over

=item :tfile => imports topen, tclose, tmkdir and trmdir

=item :temp  => imports ttempdir and ttempfile

=back

See below for a detailed description of these and other methods 
in File::Monitored.

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 DEPENDENCIES

=over

=item File::Basename

=item File::Find

=item File::Temp

=back

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package File::Monitored;

#use Moose;
use autodie qw(:all);
use File::Temp qw/ tempfile tempdir /;
use File::Basename;
use File::Find;
use Exporter;
use sigtrap 'handler' => \&_abort_cleanly, 'normal-signals', SEGV;
use strict;
use warnings;
use base qw(Exporter);

# Exported methods/tags
our @EXPORT = ();
our @EXPORT_OK = qw(topen tclose tmkdir trmdir ttempfile ttempdir);
our %EXPORT_TAGS = ('tfile' => [ qw(tclose topen tmkdir trmdir ) ],
		    'temp'  => [ qw(ttempfile ttempdir) ]);

our $VERSION = "0.2";

# Class variables
our %MONITORED_FILES = ();
our %MONITORED_DIRECTORIES = ();

=head2 topen

 Title   : topen
 Usage   : $fh = topen("data.txt")
 Function: creates a monitored file
 Returns : filehandle
 Args    : path to the file

NOTE: topen does not support the rich syntax and control arguments
of Perl's open() builtin function and always creates a new file,
erasing pre-existing files. Therefore it should never be used when
the intention is to read the contents of a file.

=cut

sub topen {
    my ($path) = @_;
    open(my $fh,"> $path") || die "Could not create monitored file $path";
    $MONITORED_FILES{$path} = 1;
    return $MONITORED_FILES{$path};
}

=head2 tmkdir

 Title   : tmkdir
 Usage   : $fh = tmkdir("data.txt")
 Function: create a monitored directory
 Returns : true if success, false if fails
 Args    : same as Perl's builtin mkdir function

=cut

sub tmkdir {
    my ($path, $mask) = @_;
    my $ret = mkdir($path, $mask);
    $MONITORED_DIRECTORIES{$path} = 1 if ($ret);
    return $ret;
}

=head2 trmdir

 Title   : trmdir
 Usage   : $fh = trmdir("data.txt")
 Function: remove a monitored directory and all its contents
 Returns : true if success, false if fails
 Args    : path to a new directory

=cut

sub trmdir {
    my ($path) = @_;
    return 0 unless (exists $MONITORED_DIRECTORIES{$path});
    finddepth(\&_remove_monitored_tree, $path);
    delete $MONITORED_DIRECTORIES{$path};
    rmdir($path);
    return 1;
}

=head2 ttempfile

 Title   : ttempfile
 Usage   : ($fh,$filename) = ttempfile()
 Function: creates a monitored temporary file
 Returns : filehandle and path to temporary file
 Args    : same as tempfile from File::Temp

See perldoc L<File::Temp>

=cut

sub ttempfile {
    my (@args) = @_;
    push(@args, 'UNLINK', 1); # Always remove temporary files
    my ($fh, $path) = tempfile(@args);
    die "Could not create temporary file" unless (defined $fh);
    $MONITORED_FILES{$path} = 1;
    return ($fh, $path);
}

=head2 ttempdir

 Title   : ttempdir
 Usage   : $fh = ttempdir()
 Function: creates a monitored temporary directory
 Returns : path to temporary directory
 Args    : none
    
NOTE: the second argument is optional but if it evaluates to
true the file is removed

=cut

sub ttempdir {
    my (@args) = @_;
    push(@args, 'CLEANUP', 1);
    my $path = tempdir(@args);
    die "Could not create temporary directory" unless (defined $path);
    $MONITORED_DIRECTORIES{$path} = 1;
    return $path;
}

=head2 INTERNAL METHODS

=head2 _abort_cleanly

 Title   : _abort_cleanly
 Usage   : _abort_cleanly()
 Function: remove all monitored files and exit
 Returns : 
 Args    : 

=cut

sub _abort_cleanly {
    my $signame = shift;

    # Removing any remaining monitored files
    foreach my $path (keys %MONITORED_FILES) {
	if (-f $MONITORED_FILES{$path}) {
	    delete $MONITORED_FILES{$path};
	    unlink($path);
	}
    }

    # Removing any remaining monitored files
    foreach my $path (keys %MONITORED_DIRECTORIES) {
	if (-d $path) {
	    finddepth(\&_remove_monitored_tree, $path);
	    delete $MONITORED_DIRECTORIES{$path};
	    rmdir($path);
	}
    }

    if ($signame) {
	die "File::Monitored caught signal $signame: $0 interrupted! Monitored files were removed...\n";
    } else {
	exit 0;
    }
}

sub _remove_monitored_tree {
    return 1 if ($_ eq '.');
    if (-d $_ && ! -l $_) {
	rmdir($_);
    } else {
	unlink($_);
    }
}

# End block to ensure clean program execution
END {
#    _abort_cleanly(0);
};

# Don't forget to return a true value to the loader!
1;
