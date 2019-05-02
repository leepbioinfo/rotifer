# POD documentation - main docs before the code

=head1 NAME

Rotifer::Tools::Rrsw::ftp - download/sync FTP site using Net::FTP 

=head1 SYNOPSIS

  # Creating a new tool script
  use Rotifer::Tools::Rrsw;
  my $rrsw = Rotifer::Tools::Rrsw->create("ftp", { name => $target });
  $rrsw->get(@ARGV);

=head1 DESCRIPTION

Rotifer::Tools::Rrsw::ftp uses Net::FTP to maintain mirrors of
ftp sites.

=head2 Net::FTP options

Options named with the prefix 'ftp_' are parsed to remove the prefix
and sent to Net::FTP. A simple example is the ftp_Passive=1 option
which activates passive mode transfers (see Net::FTP documentation).

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::Tools::Rrsw::ftp;

use strict;
use warnings;
use autodie qw(:all);
use Carp::Clan qr/^Rotifer::Tools/;
use Cwd;
use File::Basename;
use File::Spec;
use File::Copy;
use File::Listing qw(parse_dir);
use File::Path qw(make_path remove_tree);
use File::Spec;
use File::stat;
use Net::FTP;
use Moose;
use Moose::Util::TypeConstraints;
use Rotifer;
use URI;

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

=head2 options

 Usage   : Rotifer::Tools::Rrsw->new(options => { ftp_Passive => 1 })
 Function: advanced options set by rrsw or the target's configuration file
 Value   : array reference
 Builder : _default_options
 Trigger : _options_trigger

 This modules sends any option starting with 'ftp_' to the
 Net::FTP library and knows the following additional options:

=cut


subtype "HashRefWithDefaults",
     as "HashRef",
  where {
	  croak "Options attribute should be a hash reference"
	     unless (ref($_) eq "HASH");
	  if (exists $_->{"_coerced"}) {
		  delete $_->{"_coerced"};
		  return 1;
	  } else {
		  return 0;
	  }
  };

coerce "HashRefWithDefaults",
  from "HashRef",
   via {
     my $n = $_;
     my $d = &_default_options();
     map { $n->{$_} = $d->{$_} if (!exists $n->{$_}) } keys %$d;
     $n->{'_coerced'} = 1;
     return $n;
   };

has 'options' => (
    is       => "rw",
    isa      => "HashRefWithDefaults",
    lazy     => 1,
    default  => \&_default_options,
    coerce   => 1,
    );

sub _default_options {
    return { ftp_Passive => 1, follow_symlinks => 1 };
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

=cut

sub BUILD {
    my $self = shift;

    # Load configuration
    (my $package = __PACKAGE__) =~ s/Rotifer::Tools::Rrsw:://;
    my $conf = File::Spec->catfile("$Rotifer::RealScriptEtc", $self->target.".yml");
    if (exists $self->options->{'config'} && -f $self->options->{'config'}) {
	   $conf = $self->options->{config};
	   carp "Loading user-defined FTP configuration from file $conf for target ".$self->target
		  if ($self->debug > 1); 
    }
    $conf = File::Spec->catfile("$Rotifer::RealScriptEtc", "${package}.yml") unless ( -f $conf );
    my $cfg = {};
    if ( -f $conf ) {
	$cfg = Config::Any->load_files({ files           => [ $conf ],
					 use_ext         => 1,
					 flatten_to_hash => 1
				       });
	if (exists $cfg->{$conf}) {
	    $cfg = { %{$self->options}, %{$cfg->{$conf}} };
	} else {
	    die "Error while parsing $conf";
	}
    }

    # Process include/exclude lists
    foreach my $list (qw(include exclude)) {
	if (exists $cfg->{$list}) {
	    $cfg->{$list} = [ map { qr{$_} } @{$cfg->{$list}} ];
	} else {
	    $cfg->{$list} = [];
	}
    }

    # Store configuration
    $self->options($cfg);
}

=head2 execute

 Title   : execute
 Usage   : $count = $rrsw->execute("ftp://example.com/pub","/databases/src")
 Function: download/sync mirror
 Returns : list of files donwloaded
 Args    : origin and destination

=cut

sub execute {
    my ($self, $orig, $dest) = @_;
    #use Data::Dump qw(dump); print dump($self->options)."\n";

    # Prepare command-line
    if ( -d $dest ) {
	chdir($dest);
	$dest = getcwd(); # To avoid problems with symlinks
	carp "Starting script at local directory $dest" if ($self->debug);
    }

    # Connect
    my $uri = URI->new($orig);
    my %ftpopts = map { (my $n = $_) =~ s/ftp_//; (ucfirst($n),$self->options->{$_}) } grep { /^ftp_/ } keys %{$self->options};
    my $ftp = Net::FTP->new($uri->host, Debug => ($self->debug > 3 ? 1 : 0), %ftpopts);
    $ftp->hash(1,1048576) if ( -t \*STDERR ); # Show progress
    $ftp->login("anonymous",'-@anonymous');
    $ftp->binary; # Set binary mode...
    my $host = $uri->host;

    # Parse root directory
    my $remoteroot = $uri->path;
    unless ($remoteroot eq '/') {
	$remoteroot =~ s|^/+||;
	$remoteroot =~ s|/+$||;
    }

    # Go to root (local and remote)
    $ftp->cwd($remoteroot) || do { carp "Failed to cd to $remoteroot on FTP server" && return 1 };
    my $remotepwd = $ftp->pwd();
    $remotepwd =~ s|^/+|| unless ($remotepwd eq "/"); # Remove leading "/" (root)
    carp "Remote root directory set to $remotepwd (FTP root: $remoteroot)" if ($self->debug);
    my $localroot = File::Spec->catfile($dest, $host, $remoteroot);
    carp "Current local root directory is " . getcwd() if ($self->debug);

    # If $remotepwd != $remoteroot then $remoteroot is a symlink in the FTP server
    if ($remotepwd ne "/" && $remotepwd ne $remoteroot) {
	my $realroot = File::Spec->catfile($dest, $host, $remotepwd);
	if (! $self->test ) {
	    if ( -l $localroot ) {
		my $linked = Cwd::abs_path(File::Spec->rel2abs(readlink($localroot)));
		if ($linked =~ m|^$dest/$host/|) {
		    (my $rellinked = $linked) =~ s|^$dest/$host/*||;
		    if ($rellinked ne $remotepwd) {
			# Symlink must be updated: remove old linked directory and symlink
			remove_tree($linked) if ( -e $linked );
			unlink($localroot);
		    }
		} else {
		    # Link target is incompatible with rrsw!!!
		    unlink($localroot);
		}
	    } elsif ( -e $localroot ) { # $localroot exists but it is not a symlink!
		remove_tree($localroot);
	    }
	    symlink($realroot,$localroot) unless ( -l $localroot );
	}
	$localroot = $realroot;
	$remoteroot = $remotepwd;
    }

    # Prepare root of local mirror
    if ( $self->test ) {
	carp "Creating local root directory $localroot" unless ( -d $localroot );
    } else {
	make_path($localroot) unless ( -d $localroot );
	chdir($localroot);
    }
    carp "Local root directory set to $localroot" if ($self->debug);

    # Download
    my $status = $self->_process($ftp, File::Spec->catfile($dest, $host));
    $ftp->cwd("/");
    carp "Remote root directory set to " . $ftp->pwd()  if ($self->debug);
    chdir($dest); # Return to were we started
    carp "Local root directory set to " . getcwd() if ($self->debug);

    $ftp->close;
    return $status;
}

=head2 _process

 Title   : _process
 Usage   : $count = $rrsw->_process($ftp, "/blast")
 Function: evaluate and process FTP server's files/directories
 Returns : true if anything was downloaded, false otherwise
 Args    : 
  - Net::FTP object
  - root of the directory to where files should be downoaded
  - (boolean) check if target diretory matches selection rules

=cut

sub _process {
    my ($self, $ftp, $dest, $check) = @_;

    # Make sure local directory matches remote path
    my $remotepwd = $ftp->pwd;
    return 0 if ($check && $self->_ignore($remotepwd));
    (my $localpwd = $dest . $remotepwd) =~ s|//|/|g;

    # Process the contents of the current directory
    my $status = 0;  # Was anything downloaded?
    foreach my $file (parse_dir($ftp->dir)) { # Note that 'dir' ignores hidden files ("\.*")!
	my ($name, $type, $size, $mtime, $mode) = @$file; # Remote file/directory properties
	next if ($name eq '.' || $name eq '..'); # Ignore special files

	# Store current remote target
	(my $remotepath = "$remotepwd/$name") =~ s|//|/|g;
	carp "Current target is $remotepath" if ($self->debug > 1);
        next if $self->_ignore($remotepath);

	# Store full URL
	(my $uri = "ftp://" . $ftp->host . $remotepath) =~ s|//|/|g;

	# Scan remote directories that were not explicitly excluded ('next if' statement above)
	if ($type eq 'd') {
	    my $ok = $ftp->cwd($name);
	    if (!$ok) {
	        carp "Could not access diretory $name at ".$remotepath if ($self->debug > 1);
	        next;
            }
	    carp "Remote directory set to " . $ftp->pwd() if ($self->debug > 2);
	    my $localstatus = $self->_process($ftp, $dest, 1);
	    if ($localstatus && $remotepath ne $ftp->pwd) {
		(my $orig = "$dest/" . $ftp->pwd) =~ s|//|/|g;
		carp "Adding symlink: $localpwd/$name -> $orig" if ($self->debug > 1);
		$self->_mkcwd($localpwd);
		if (!$self->test) {
		    remove_tree($name) if ( -e $name || -l $name );
		    symlink($orig, $name);
		}
	    }
	    $ftp->cwd($remotepwd);
	    carp "Remote directory set to " . $ftp->pwd() if ($self->debug > 2);
	    $status ||= $localstatus;
	}

        # Process symlinks to files
	elsif ($type =~ /^l +(\S+)/) {
            my $linked = $1;
	    #$links{"$localpwd/$name"} = $linked;
            my @contents = grep { !/^\.+$/ } $ftp->ls($name);
	    my $index = index($contents[0],"/") if (scalar @contents);
	    my $localstatus = 0;
	    if (defined $index && $index != -1 && $name eq substr($contents[0],0,$index)) {
		carp "$remotepath is a symlink to directory $linked" if ($self->debug > 1);
		my $ok = $ftp->cwd($linked);
		if (!$ok) {
		    carp "Could not access diretory $linked at " if ($self->debug > 1);
		    next;
		}
		carp "Remote directory set to " . $ftp->pwd() if ($self->debug > 2);
		$localstatus = $self->_process($ftp, $dest, 1);
		$ftp->cwd($remotepwd);
		carp "Remote directory set to " . $ftp->pwd() if ($self->debug > 2);
	    } else {
		    $localstatus = 1;
	    }
	    if ($localstatus) {
		    (my $orig = "$dest/$linked") =~ s|//|/|g;
		    carp "Adding symlink: $localpwd/$name -> $orig" if ($self->debug > 1);
		    $self->_mkcwd($localpwd);
		    if (!$self->test) {
			remove_tree($name) if ( -e $name || -l $name );
			symlink($orig, $name);
			$status ||= $localstatus;
		    }
	    }
	}

	# Regular files
	elsif ($type eq "f") {
	    $self->_mkcwd($localpwd);

	    # Evaluate file properties
	    my $filestate = 0;
            my $stat = stat($name);
	    if (defined $stat) {
	        $filestate += 1 if ($stat->size != $ftp->size($name));
		$filestate += 2 if ($stat->mtime < $ftp->mdtm($name));
		next if (!$filestate); # $filestate == 0
	     	if ($self->debug > 1) {
                    if ($filestate == 1) {
			carp "Size change detected for file $uri";
		    } elsif ($filestate == 2) {
			carp "Newer upstream file $uri";
		    } elsif ($filestate == 3) {
			carp "Size change detected for newer upstream file $uri"
		    }
		}
            }

	    # Download or report new/changed
	    if ($self->test) {
	        carp "get ($filestate) $uri";
		$status = 1;
	    } else {
		carp "retrieving ($filestate) $uri";
		my $localpath = $ftp->get($name);
	     	carp "downloaded copy of $uri to $localpath" if ($self->debug > 2);
	    	$status ||= defined $localpath;
	    }
	} # elsif ($type eq "f")

	# File is not a directory, symlink or regular file: what????
	else {
	    carp "Ignoring unknown file type: $remotepath ($type)";
	}
    } # foreach my $file (parse_dir($ftp->dir))

    return $status;
}

=head2 _mkcwd

 Title   : _mkcwd
 Usage   : $count = $rrsw->_activate_diretory($path)
 Function: create and change to a local directory
 Returns : 

=cut

sub _mkcwd {
    my ($self, $path) = @_;
    return if ($path eq getcwd());
    make_path($path) unless ( -d $path || $self->test);
    carp "Changing local directory from ".getcwd()." to $path" if ($self->debug > 2);
    if ( -d $path ) {
        chdir($path);
        carp "Local directory set to " . getcwd() if ($self->debug > 2);
    } elsif (!$self->test) {
        croak "missing destination directory $path";
    }
}

=head2 _ignore

 Title   : _ignore
 Usage   : $count = $rrsw->_ignore($path)
 Function: test whether $path matches regular 
           expression inclusion/exclusion patterns
 Returns : 

=cut

sub _ignore {
	my ($self, $path) = @_;

	# Check if this target matches any inclusion pattern
	my @inclMatch = grep { $path =~ m{$_} } @{$self->options->{include}};
        if (!scalar @inclMatch) {
            carp "Target $path doesn't match any filters (i.e. include regexps)." if ($self->debug > 2);
	    return 1;
        } else {
            carp "Target $path matches the following filters: @inclMatch" if ($self->debug > 2);
	}

	# Ignore excluded files/diretories
	my @exclMatch  = grep { $path =~ m{$_} } @{$self->options->{exclude}};
	if (scalar @exclMatch) {
		carp "Target $path matches the following excluded patterns: @exclMatch" if ($self->debug > 2);
		return 1;
	}

	# Accept
	return 0;
}

__PACKAGE__->meta->make_immutable;
1;
