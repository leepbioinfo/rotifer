# POD documentation - main docs before the code

=head1 NAME

Rotifer::Tools::Rrsw::base - base class for rrsw plugin modules

=head1 SYNOPSIS

  # Creating a new tool script
  use Rotifer::Tools::Rrsw;
  my $rrsw = Rotifer::Tools::Rrsw->create("base");
  $rrsw->execute($indir, $outdir, @args);

=head1 DESCRIPTION

This module provides basic methods to processes the contents of
a directory tree and use the output files generated to add/replace
the contents of another directory.

Plugins should inherit from this class and overload the method
_execute(). subroutines to each file within the input directory
and saving the output to a similar directory tree in another
directory.

Once _execute() is run by calling $rrsw->execute(), files generated
by the child class will be moved to outiput directory. Older 
versions of these files, if found in $outdir, will be atomically
replaced by new the ones.

The list of subroutines given as argument 

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::Tools::Rrsw::base;

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
use Moo;
use Rotifer;
use Scalar::Util qw(blessed);

=head2 target

 Usage   : Rotifer::Tools::Rrsw::base->new(target => $target)
 Function: basename of the target
 Value   : string
 Required: 1

=cut

has 'output' => (
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

 Usage   : Rotifer::Tools::Rrsw::base->new(options => { gzip => "/usr/bin/gzip" })
 Function: module configuration
 Value   : hash reference as returned by Config::General

=cut

has 'options' => (
    is       => "rw", 
    isa      => "HashRef",
    );

=head2 callback

 Usage   : Rotifer::Tools::Rrsw::base->new(callback => [{ coderef => $s, args => \@a }])
 Function: list of subroutines to apply to each target
 Value   : hash reference as returned by Config::General

=cut

has 'callback' => (
    is       => "rw", 
    isa      => "ArrayRef[HashRef]",
    default  => sub { [] },
    );

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

=head2 execute

 Title   : execute
 Usage   : $count = $rrsw->execute($indir, $dest)
 Function: callback downloaded files
 Returns : list of files eachnew
 Args    : input file(s)/directory and output directory

=cut

sub execute {
    my ($self, $indir, $outdir) = @_;

    # Prepare to work in a temporary directory
    my $tmpdir;
    unless ($self->test) {
	mkdir($outdir) if ( ! -e $outdir );
	chdir($outdir) || croak "Target directory $outdir is unaccessible";
	$tmpdir = tempdir(".".$self->target.".XXXXXXXX");
	chdir($tmpdir) || croak "Unable to access temporary processed directory $tmpdir";
	map { mkdir($_) } qw(new old);
	chdir("new")|| croak "Unable to access temporary processed directory $tmpdir/new";
    }

    # Apply processing subroutines to each source path
    my @subroutines = @{$self->callback};
    for (my $i=0; $i<=$#subroutines; $i++) {
       if (!exists $subroutines[$i]->{'coderef'}) {
          if (exists $subroutines[$i]->{'code'}) {
             $subroutines[$i]->{'coderef'} = eval $subroutines[$i]->{'code'};
             croak "Error while parsing rule #$i:\n$@\nCode:\n".$subroutines[$i]->{'code'}."\n" if ($@);
	  } else {
            croak "Incomplete callback subroutine definition for target ".$self->target;
	  }
       }
       my @args = $subroutines[$i]->{'args'};
       @args = ref($args[0]) eq "HASH" ? %{$args[0]} : @{$args[0]} if (defined(ref($args[0])));
       my $status = $subroutines[$i]->{'coderef'}->($self,$indir,$outdir,@args);
    }

    # Move new files to $outdir almost atomically.
    # TODO: make it atomic (File::Transaction::Atomic?)
    unless ($self->test) {
	foreach my $patt ("*",".*") {
	    foreach my $eachnew (grep { !/^\.+$/ } glob($patt)) {
		if ( -e "$outdir/$eachnew" ) {
		    my $success = rename("$outdir/$eachnew","../old/$eachnew");
		    croak "can't replace file/directory $outdir/$eachnew" if (!$success);
		}
		my $success = rename($eachnew,"$outdir/$eachnew");
		croak "can't replace file/directory $outdir/$eachnew" if (!$success);
	    }
	}
	chdir($outdir);
	remove_tree($tmpdir);
    }

    chdir($Rotifer::IWD);
    return;
}

__PACKAGE__->meta->make_immutable;
1;
