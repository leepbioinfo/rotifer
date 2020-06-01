# POD documentation - main docs before the code

=head1 NAME

Rotifer::Tools::Rrsw::unpack_blast - unpack and format Blast databases

=head1 SYNOPSIS

  # Creating a new tool script
  use Rotifer::Tools::Rrsw;
  my $rrsw = Rotifer::Tools::Rrsw->create("unpack_blast", { name => $target });
  $rrsw->get(@ARGV);

=head1 DESCRIPTION

Rotifer::Tools::Rrsw::unpack_blast unpacks Blast database files and prepares
these files to be used in searches by BLAST and HMMER 3 (using blastdbcmd
and esl-sfetch).

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::Tools::Rrsw::unpack_blast;

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

 Usage   : Rotifer::Tools::Rrsw::unpack_blast->new(target => $target)
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

 Usage   : Rotifer::Tools::Rrsw::unpack_blast->new(options => { gzip => "/usr/bin/gzip" })
 Function: module configuration
 Value   : hash reference as returned by Config::General
 Builder : _default_options

=cut

# Inherited from Rotifer::Tools::Rrsw::unpack

=head2 executable

 Usage   : Rotifer::Tools::Rrsw::unpack_blast->new(executable => { gzip => "/usr/bin/gzip" })
 Function: path to tar, gzip and other executable
 Value   : hash reference with program names as keys
 Builder : _default_executable (uses File::Which)

=cut

sub _default_executable {
    my $self = shift;
    my $paths = $self->next::method();
    my $path = which('blastdbcmd');
    $paths->{blastdbcmd} = $path if (defined $path);
    $path = which('esl-sfetch');
    $paths->{"esl-sfetch"} = $path if (defined $path);
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

    # Detect appropriate executables
    my $blastdbcmd = $self->executable->{'blastdbcmd'};
    die "Could not locate blastdbcmd. You need to install the NCBI's C++ Toolkit or BLAST++."
	unless (defined $blastdbcmd);
    my $eslsfetch = $self->executable->{'esl-sfetch'};
    die "Could not locate esl-sfetch. You should install HMMER with the configure flag --with-easel."
	unless (defined $eslsfetch);

    # Unpack first
    $self->next::method($source, $unpack);

    # Format unpacked databases
    mkdir($unpack) if ( ! -d $unpack );
    chdir($unpack);

    # Process options
    my @exclude = @{$self->options->{'exclude'}} if (exists $self->options->{'exclude'});

    # Identify databases
    my $dbs;
    run([ $blastdbcmd, "-recursive", "-remove_redundant_dbs", "-list", ".", "-list_outfmt", '%f', '-dbtype', 'prot' ], \undef, \$dbs);

    # Run
    my $status = 0;
    foreach my $db (split(/\n/,$dbs)) {
	next if ($db =~ /\.\d\d$/); # Avoid fragments
	next if (grep { $db =~ /$_/ } @exclude);
	my @cmd1 = ([ $blastdbcmd, "-db", $db, "-entry", "all", "-out", "${db}.fa" ]);
	my @cmd2 = ([ "perl", "-i", "-pe", 's/\-/X/go if (substr($_,0,1) ne ">")', "${db}.fa" ]);
	my @cmd3 = ([ $eslsfetch, "--index", "${db}.fa" ]);
	if ($self->test) {
	    print join(" ",@{$cmd1[0]}),"\n";
	    print join(" ",@{$cmd2[0]}),"\n";
	    print join(" ",@{$cmd3[0]}),"\n";
	} else {
	    $status = run(@cmd1);
	    $status = run(@cmd2);
	    $status = run(@cmd3) if ( -f "${db}.fa" );
	}
    }

    chdir($Rotifer::IWD);
    return $status;
}

__PACKAGE__->meta->make_immutable;
1;
