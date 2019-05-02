# POD documentation - main docs before the code

=head1 NAME

Rotifer::Tools::Rrsw::dump_blast - unpack and format PFam databases

=head1 SYNOPSIS

  # Creating a new tool script
  use Rotifer::Tools::Rrsw;
  my $rrsw = Rotifer::Tools::Rrsw->create("dump_blast", { name => $target });
  $rrsw->get(@ARGV);

=head1 DESCRIPTION

Rotifer::Tools::Rrsw::dump_blast unpacks BLAST database files and prepares
these files to be used in searches by HMMER 3 (using hmmpress).

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::Tools::Rrsw::dump_blast;

use strict;
use warnings;
use autodie qw(:all);
use Carp::Clan qr/^Rotifer::Tools::Rrsw/;
use File::Basename;
use File::Find;
use File::Which;
use IPC::Run qw(run);
use Moose;
use Rotifer;
extends 'Rotifer::Tools::Rrsw::unpack';

=head2 target

 Usage   : Rotifer::Tools::Rrsw::dump_blast->new(target => $target)
 Function: basename of the target
 Value   : string
 Required: 1

=cut

# Inherited from Rotifer::Tools::Rrsw::unpack

=head2 test

 Usage   : Rotifer::Tools::Rrsw->new(test => 1)
 Function: run in test mode
 Value   : Boolean

=cut

# Inherited from Rotifer::Tools::Rrsw::unpack

=head2 options

 Usage   : Rotifer::Tools::Rrsw::dump_blast->new(options => { gzip => "/usr/bin/gzip" })
 Function: module configuration
 Value   : hash reference as returned by Config::General
 Builder : _default_options

=cut

# Inherited from Rotifer::Tools::Rrsw::unpack

=head2 executable

 Usage   : Rotifer::Tools::Rrsw::dump_blast->new(executable => { gzip => "/usr/bin/gzip" })
 Function: path to tar, gzip and other executable
 Value   : hash reference with program names as keys
 Builder : _default_executable (uses File::Which)

=cut

sub _default_executable {
    my $self = shift;
    my $paths = $self->next::method();
    my $path = which('blastdbcmd');
    $paths->{blastdbcmd} = $path if (defined $path);
    return $paths;
}

=head2 METHODS

=head2 execute

 Title   : execute
 Usage   : $count = $rrsw->execute($source, $dest)
 Function: unpack downloaded files
 Returns : list of files unpacked
 Args    : input file(s)/directory and output directory

=cut

sub execute {
    my ($self, $source, $unpack) = @_;

    # Prepare to work in a temporary directory
    my $tmpdir;
    unless ($self->test) {
	mkdir($unpack) if ( ! -e $unpack );
	chdir($unpack);
	$tmpdir = tempdir(".".$self->target.".XXXXXXXX");
	chdir($tmpdir);
	map { mkdir($_) } qw(new old);
	chdir("new");
    }

    # Unpack everything here
    my @source = ($source);
    if ( -d $source ) {
	find({ no_chdir => 1,
	       wanted   => sub { push(@source, $File::Find::name) },
	     }, $source);
    }
    map { $self->unpack_each($source, $_) } @source;

    # Move new files to $unpack almost atomically. TODO: make it atomic (File::Transaction::Atomic?)
    unless ($self->test) {
	foreach my $patt ("*",".*") {
	    foreach my $unpacked (grep { !/^\.+$/ } glob($patt)) {
		if ( -e "$unpack/$unpacked" ) {
		    my $success = rename("$unpack/$unpacked","../old/$unpacked");
		    croak "can't replace file/directory $unpack/$unpacked" if (!$success);
		}
		my $success = rename($unpacked,"$unpack/$unpacked");
		croak "can't replace file/directory $unpack/$unpacked" if (!$success);
	    }
	}
	chdir($unpack);
	remove_tree($tmpdir);
    }

    chdir($Rotifer::IWD);

    # Format unpacked databases
    mkdir($unpack) if ( ! -d $unpack );
    chdir($unpack);

    # Detect appropriate executable
    my $blastdbcmd = $self->executable->{'blastdbcmd'};
    die "Could not locate blastdbcmd. You need to install the NCBI's C++ Toolkit or BLAST++."
	unless (defined $blastdbcmd);

    # Identify databases
    my $dbs;
    run([ $blastdbcmd, "-recursive", "-remove_redundant_dbs", "-list", ".", "-list_outfmt", '%f' ], \undef, \$dbs);

    # Run
    my $status = 0;
    foreach my $db (split(/\n/,$dbs)) {
	next if ($db =~ /\.\d\d$/); # Avoid fragments
	my @cmd = ([ $blastdbcmd, "-db", $db, "-entry", "all", "-out", "${db}.fa" ]);
	if ($self->test) {
	    print join(" ",@{$cmd[0]}),"\n";
	} else {
	    $status = run(@cmd);
	}
    }

    chdir($Rotifer::IWD);
    return $status;
}

__PACKAGE__->meta->make_immutable;
1;
