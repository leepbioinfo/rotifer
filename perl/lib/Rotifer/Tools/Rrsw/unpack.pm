# POD documentation - main docs before the code

=head1 NAME

Rotifer::Tools::Rrsw::unpack - simple unpacker for multiple file formats 

=head1 SYNOPSIS

  # Creating a new tool script
  use Rotifer::Tools::Rrsw;
  my $rrsw = Rotifer::Tools::Rrsw->create("unpack", { name => $target });
  $rrsw->get(@ARGV);

=head1 DESCRIPTION

Rotifer::Tools::Rrsw::unpack unpacks files in several different formats
like Tar, Gzip, etc.

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::Tools::Rrsw::unpack;

use strict;
use warnings;
use autodie qw(:all);
use Carp::Clan qr/^Rotifer::Tools::Rrsw/;
use Config::Any;
use File::Basename;
use File::Find;
use File::Path qw(remove_tree make_path);
use File::Spec;
use File::Temp qw(tempdir);
use File::Which;
use IPC::Run qw(run);
use Moose;
use Rotifer;
use Scalar::Util qw(blessed);

=head2 target

 Usage   : Rotifer::Tools::Rrsw::unpack->new(target => $target)
 Function: basename of the target
 Value   : string
 Required: 1

=cut

has 'target' => (
    is       => "rw", 
    isa      => "Str",
    required => 1,
    );

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

=head2 options

 Usage   : Rotifer::Tools::Rrsw::unpack->new(options => { gzip => "/usr/bin/gzip" })
 Function: module configuration
 Value   : hash reference as returned by Config::General
 Builder : _default_options

=cut

has 'options' => (
    is       => "rw", 
    isa      => "HashRef",
    lazy     => 1,
    builder  => '_default_options',
    );

sub _default_options {
    return {
	"tar" => {
	    command_line_options =>  [ "-xv" ],
	},
	"gzip" => {
	    command_line_options =>  [ "-dc" ],
	},
	"bzip2" => {
	    command_line_options =>  [ "-dc" ],
	},
    };
}

=head2 executable

 Usage   : Rotifer::Tools::Rrsw::unpack->new(executable => { gzip => "/usr/bin/gzip" })
 Function: path to tar, gzip and other executable
 Value   : hash reference with program names as keys
 Builder : _default_executable (uses File::Which)

=cut

has 'executable' => (
    is       => "rw", 
    isa      => "HashRef[Str]",
    lazy     => 1,
    builder  => '_default_executable',
    );

sub _default_executable {
    my %paths = ();
    foreach my $exe (qw(tar gzip bzip2 unzip)) {
	my $path = which($exe);
	$paths{$exe} = $path if (defined $path);
    }
    return \%paths;
}

=head2 suffixes

 Usage   : Rotifer::Tools::Rrsw::unpack->new(suffixes => { '.tar.gz' => "tar" })
 Function: map of suffixes to unpacking programs
 Value   : hash reference with filename suffixes as keys
 Builder : _default_executable (uses File::Which)

=cut

has 'suffixes' => (
    is       => "rw", 
    isa      => "HashRef[Str]",
    lazy     => 1,
    builder  => '_default_suffixes',
    );

sub _default_suffixes {
    my ($self) = @_;
    my $ext = {
	'.tar.gz' => "tar",
	'.tgz'    => "tar",
	'.tar'    => "tar",
	'.gz'     => "gzip",
	'.zip'    => "unzip",
	'.bz2'    => "bzip2",
    };
    return $ext;
}

=head2 METHODS

=cut

sub BUILD {
    my $self = shift;

    # Load configuration
    (my $PACKAGE = __PACKAGE__) =~ s/Rotifer::Tools::Rrsw:://;
    foreach my $target ($PACKAGE, $self->target) {
	my $conf = File::Spec->catfile("$Rotifer::RealScriptEtc","${target}.yml");
	next unless ( -f $conf );
	my $cfg = Config::Any->load_files({ files           => [ $conf ],
					    use_ext         => 1,
					    flatten_to_hash => 1
					  });
	if ($cfg->{$conf}) {
	    $cfg = { %{$self->options}, %{$cfg->{$conf}} };
	    $self->options($cfg);
	}
    }

    #use Data::Dump qw(dump); die dump($self->options);
}

sub _compile {
  my @FILTERS = ();
  foreach (@_) {
      my $userdir = Rotifer->userdir;
      my $sysdir = Rotifer->sysdir;
      my $ref = $_;
      my $rule = $_;
      if ( -f "$userdir/lib/rrsw/$rule" ) {
          $rule = "$userdir/lib/rrsw/$rule";
      } elsif ( -f "$sysdir/lib/rotifer/rrsw/$rule" ) {
          $rule = "$sysdir/lib/rotifer/rrsw/$rule";
      }
      if ($rule =~ /.sh$/) { # Shell script!
          $ref = $rule;
      } else {
          if (-f "$rule") {
              open(RULE,"<$rule") || die "Could not load rules from file $rule";
              $ref = join(" ",map { chomp; $_ } grep { !/^\#/ } <RULE>);
              close(RULE);
          }
	  $ref = eval "$ref";
          die "Error while compiling output filter $rule\n$@" if ($@);
      } 
      push(@FILTERS, $ref);
  }
  return @FILTERS
}

sub _run_shell_script {
  my ($self, $source, $unpack, $command) = @_;
  use IPC::Run qw(run);
  my $exitcode = run([ $command, $self->test, $source, $unpack ]);
  $exitcode = $exitcode >> 8;
  return $exitcode;
}

=head2 execute

 Title   : execute
 Usage   : $count = $rrsw->execute($source, $dest)
 Function: unpack downloaded files
 Returns : list of files unpacked
 Args    : input file(s)/directory and output directory

=cut

sub execute {
    my ($self, $source, $unpack) = @_;

    # Process callback processsors
    my @callback = ();
    if (exists $self->options->{'callback'}) {
        @callback = _compile(@{$self->options->{'callback'}});
    }

    # Prepare to work in a temporary directory
    my $tmpdir;
    if ($self->test) {
        chdir($unpack)
    } else {
	mkdir($unpack) if ( ! -e $unpack );
	chdir(dirname($unpack));
	$tmpdir = tempdir(".".$self->target.".XXXXXXXX");
	chdir($tmpdir);
	#map { mkdir($_) } qw(new old);
	#chdir("new");
    }

    # Unpack everything here
    my @source = ($source);
    if ( -d $source ) {
	find({ no_chdir => 1,
	       wanted   => sub { push(@source, $File::Find::name) },
	     }, $source);
    }
    map { $self->unpack_each($source, $_) } @source;

    # Call callback function to process files
    if (scalar @callback) {
	# Running Perl callback functions
        my @code = grep { ref($_) eq "CODE" } @callback;
	if (scalar @code) {
            find({ no_chdir => 1,
                   wanted   => sub {
                      if ($File::Find::name ne ".") {
                          foreach my $method (@code) {
                             $method->($self, $source, $unpack);
                          }
                      }
                  }
              }, $self->test ? $unpack : ".");
        }

	# Running external shell scripts
	my $i = 0;
	my @scripts = grep { $_ =~ /\.sh$/ } @callback;
	foreach my $script (@scripts) {
	   _run_shell_script($self, $source, $unpack, $script, $i);
	   $i++;
        }
    }

    # Move new files to $unpack almost atomically. TODO: make it atomic (File::Transaction::Atomic?)
    unless ($self->test) {
    #	foreach my $patt ("*",".*") {
    #	    foreach my $unpacked (grep { !/^\.+$/ } glob($patt)) {
    #		if ( -e "$unpack/$unpacked" ) {
    #		    my $success = rename("$unpack/$unpacked","../old/$unpacked");
    #		    croak "can't replace file/directory $unpack/$unpacked" if (!$success);
    #		}
    #		my $success = rename($unpacked,"$unpack/$unpacked");
    #		croak "can't replace file/directory $unpack/$unpacked" if (!$success);
    #	    }
    #	}
    #	remove_tree($tmpdir);
        use POSIX 'strftime';
        my $date = strftime '%Y%m%d_%H%M', localtime;
        chdir(dirname($unpack));
	sleep 60 if ( -d "${unpack}.$date" );
	chmod(0775, $tmpdir)
        rename($unpack,"${unpack}.$date");
        rename($tmpdir,$unpack);
    }

    chdir($Rotifer::IWD);
    return;
}

=head2 unpack_each

 Title   : unpack_each
 Usage   : $count = $rrsw->unpack_each($file)
 Function: unpack a file using tar, gzip or bzip2
 Returns : exit status of backend program
 Args    : input file path

=cut

sub unpack_each {
    my ($self,$source,$file) = @_;
    return unless ( -f $file || -l $file );
    my $status;

    # Detect appropriate executable
    my ($cmd,$ext) = $self->get_executable_by_extension($file);
    my $basename = basename($file, $ext || "");

    # For all methods but tar
    if (!defined $cmd || $cmd ne "tar") { # Why do I ignore tar?
	my $dir = dirname($file);
	unless ("$dir/" eq $source || $dir eq $source) {
	    $dir =~ s|^${source}/?||;
	    make_path($dir) unless (-d $dir || $self->test);
	    $basename = "$dir/$basename";
	}
    }

    # Prepare command line
    if (defined $cmd) {
	# Load options
	my @opts = ();
	if (exists $self->options->{$cmd} && exists $self->options->{$cmd}{command_line_options}) {
	    @opts = @{ $self->options->{$cmd}{command_line_options} };
	}

	# Process with tar
	my @cmd = (defined $cmd && exists $self->executable->{$cmd} ? $self->executable->{$cmd} : $cmd);
	if ($cmd eq "tar") {
	    @cmd = ([ @cmd, "-f", $file, @opts ]);
	}

	# Process with gzip or bzip2
	elsif ($cmd eq "gzip" || $cmd eq "bzip2") {
	    @cmd = ([ @cmd, $file, @opts ], \undef, ">$basename");
	}

	# Run command for each source
	if ($self->test) {
	    print join(" ",map { ref $_ eq 'ARRAY' ? @$_ : ref $_ eq 'SCALAR' ? $$_ || "" : $_ } @cmd),"\n";
	} else {
	    run(@cmd);
	}
    } # if (defined $cmd)

    # If you don't know what to do, symlink!
    else { # elsif ($cmd eq "link" || $cmd eq "ln") {
	$basename .= $ext if (defined $ext); # Restore simplest basename
	if ($self->test) {
	    print "ln -s $file $basename\n"; # Just show for the user what would be done
	} else {
	    symlink($file, $basename) unless (-e $basename);
	}
    }

    return $?;
}

=head2 get_executable_by_extension

 Title   : get_executable_by_extension
 Usage   : $count = $rrsw->get_executable_by_extension($file)
 Function: identify appropriate program for each file based
           on the table of known filename suffixes
 Returns : executable name and suffix (two strings)
 Args    : filename

=cut

sub get_executable_by_extension {
    my ($self, $filename) = @_;
    my ($cmd,$ext);
    foreach my $suffix (sort { length($b) <=> length($a) } keys %{$self->suffixes}) {
	if ($filename =~ /$suffix$/) {
	    $cmd = $self->suffixes->{$suffix};
	    $ext = $suffix;
	    last;
	}
    }
    return ($cmd, $ext);
}

__PACKAGE__->meta->make_immutable;
1;
