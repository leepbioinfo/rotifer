# POD documentation - main docs before the code

=head1 NAME

Rotifer::DBIC::AnnotationDB::Parser::ontology - 
 a bioperl-based parser for ontologies

=head1 SYNOPSIS

  # Using this parser with rotiferDB command line app

  rotiferDB -if ontology te.obo

  # Creating a new parser
  use Rotifer::DBIC::Parser;
  my $parser = Rotifer::DBIC::Parser->create("ontology");
  $parser->load(@ARGV);

=head1 DESCRIPTION

Rotifer::DBIC::Parser::ontology will load ontologies using Bioperl-db.
Work is on the way to remove BioSQL's dependencies.

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 DEPENDENCIES

=over

=item Moose

=item Rotifer::DBIC

=item Carp::Clan

=back

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::DBIC::AnnotationDB::Parser::ontology;

use strict;
use warnings;
use autodie qw(:all);
use Carp::Clan qr/^Rotifer::DBIC::AnnotationDB/;
use File::Which;
use IPC::Run qw(run);
use Moose;
with 'Rotifer::DBIC::AnnotationDB::Role::ParserRole';

=head2 ATTRIBUTES / ACCESSORS

This section list attributes not defined by the fundamental AnnotationDB parser's role.

=head2 isa_transaction_manager

 Usage   : $writer->isa_transaction_manager(1)
 Function: true if the parser will control all of its transactions
 Value   : boolean
 Default : true

=cut

has '+isa_transaction_manager' => (default => 1 );

=head2  options

 Usage   : Rotifer::DBIC::AnnotationDB::Parser->new(options => { '-format' => "obo" })
 Function: bp_load_ontology.pl commend lline options
 Value   : array reference
 Default : [ --lookup --updobsolete --safe --computetc --format obo ]

=cut

has 'options'   => (
    is       => "rw",
    isa      => "ArrayRef",
    default  => sub { [qw/ --lookup --updobsolete --safe --computetc --format obo /] },
    );

=head2  executable

 Usage   : Rotifer::DBIC::AnnotationDB::Parser->new(executable => "bp_load_ontology.pl")
 Function: path or name of the bp_load_ontology.pl script that is bundled to BioSQL

           If only the name is given an attempt is made to locate the script using
           File::Which

 Value   : string
 Default : bp_load_ontology.pl

=cut

has 'executable'   => (
    is       => "rw",
    isa      => "Str",
    builder  => "_default_executable",
    trigger  => \&_executable_trigger,
    );

sub _default_executable {
    my $loader = which("bp_load_ontology.pl");
    die "Could not find the bp_load_ontology.pl" unless (defined $loader);
    return $loader;
}

sub _executable_trigger {
    my ($self, $new, $old) = @_;
    if (defined $new) {
	my $path = which($new);
	die "Unable to locate script $new" unless (defined $path || -f $new);
    }
}

=head2 load

 Title   : load
 Usage   : $full_path = $parser->load("data.txt")
 Function: process and load data from all input files
 Returns : 
 Args    : list of file names

=cut

sub load {
    my ($self, @input) = @_;

    # Ugly hack for BioSQL compatibility on PostgreSQL
    #  until I move away from BioSQL
    my %rules = ();
    if ($self->schema->storage->sqlt_type eq 'PostgreSQL') {
	$self->schema->storage->debug($self->schema->debug);
	my $guard = $self->schema->txn_scope_guard;
	%rules = $self->schema->storage->dbh_do(
	    sub {
		my ($storage, $dbh, @args) = @_;

		my %fixed = ();
		my %biosql = map { ($_,1) } qw(dbxref dbxref_qualifier_value ontology term term_dbxref term_path term_relationship term_synonym);
		foreach my $source_name ($self->schema->sources) {
		    my $source = $self->schema->source($source_name);
		    my $name = $source->name;
		    next if (!exists $biosql{$name});
		    next if (exists $fixed{$name}); # Why does this loop process a source more than once?

		    # Add rules to ignore attempts to violate unique constraints of BioSQL tables
		    my %unique  = $source->unique_constraints;
		    my $nofKeys = scalar(keys %unique);
		    my $i = 1;
		    my $pkey = $unique{primary}->[0];
		    foreach my $c (keys %unique) {
			next if ($nofKeys > 1 && $c eq "primary");
			$pkey = $unique{$c}->[0] unless (defined $pkey);
			my $cond = "SELECT $pkey FROM ${name} WHERE " . join(" AND ",map { "$_ = new.$_" } @{$unique{$c}});
			my $q = "CREATE RULE rule_${name}_i$i AS ON INSERT TO ${name} WHERE ($cond) IS NOT NULL DO INSTEAD NOTHING";
			print "$q\n" if $self->schema->debug;
			$dbh->do($q);
			$rules{"rule_${name}_i$i"} = $name;
			$i++;
		    }

		    $fixed{$name} = 1;
		}

		return %rules;
	    });
	$guard->commit;
    } # if ($self->schema->storage->sqlt_type eq 'PostgreSQL')

    # Load
    my @connect_info = @{ $self->schema->storage->connect_info };
    my @cmd = ($self->executable(), "--dsn", $connect_info[0], @{ $self->options });
    push(@cmd, "--dbuser", $connect_info[1]) if (defined $connect_info[1]);
    push(@cmd, "--dbpass", $connect_info[2]) if (defined $connect_info[2]);
    push(@cmd, @input);
    print join(" ",@cmd),"\n" if ($self->schema->debug);
    my $stat = run(@cmd);

    # Drop rules for BioSQL compatibility in PostgreSQL
    my $guard = $self->schema->txn_scope_guard;
    if (scalar(keys %rules)) {
	my @stuff = $self->schema->storage->dbh_do(
	    sub {
		my ($storage, $dbh, @args) = @_;
		foreach my $rule (keys %rules) {
		    my $q = "DROP RULE IF EXISTS $rule ON $rules{$rule}";
		    print "$q\n" if $self->schema->debug;
		    $dbh->do($q);
		}
	    });
    }

    $self->schema->storage->dbh_do(sub { my ($s, $dbh, @a) = @_; $dbh->do("UPDATE term SET is_obsolete = '0' WHERE is_obsolete IS NULL") });
    $guard->commit;
    return $stat;
}

__PACKAGE__->meta->make_immutable;
1;
