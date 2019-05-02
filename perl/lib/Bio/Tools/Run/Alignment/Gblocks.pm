# $Id: Gblocks.pm,v 1.0 2006/05/04 18:02:00 rfsouza Exp $
#
# BioPerl module for Bio::Tools::Run::Align::Gblocks
#
# Cared for by Robson Francisco de Souza <rfsouza-AT-cecm_DOT_usp_DOT_br>
#
# Copyright Robson Francisco de Souza
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Run::Align::Gblocks - Wrapper aroud the Gblocks program

=head1 SYNOPSIS

  # Laoding the module
  use Bio::Tools::Run::Align::Gblocks;


  # There are basically two ways to use this wrapper.


  # The simplest one: build wrapper, name an alignment file
  # and run!
  my $wrapper = Bio::Tools::Run::Align::Gblocks->new(-alignment_file=>'a.phy');
  my ($rc, $hash) = $wrapper->run();


  # Or, a little more complicated: 
  # 1) load an alignment
  # 2) pass some arguments to the wrapper and run


  # Load an alignment
  use Bio::AlignIO;
  my $aio = Bio::AlignIO->new(-file=>'input.aln',
                              -format=>'clustalw');
  my $align = $aio->next_aln;

  # Call the wrapper to analyze the alignment
  my $wrapper = new Bio::Tools::Run::Align::Gblocks(-alignment=>$align);

  # Choose to estimate the alpha parameter for among site 
  # rate variation and make no bootstrap replicates
  $wrapper->set_parameters();

  # Running Gblocks
  my ($rc, $hash) = $wrapper->run();


  # Whatever the way you pass an alignment to the wrapper
  # you need to print something from its output.

  # Printing the best ML tree found:
  use Bio::AlignIO;
  my $tio = Bio::AlignIO->new(-file=>'>output.aln',
                              -format=>'phylip');
  $tio->write_aln($hash->{alignment});


=head1 DESCRIPTION

This module implements a wrapper around the Gblocks program
for selecting conserved clumns in multiple sequence alignments.
It accepts both alignment files and alignment objects as input,
through the methods alignment_file() and alignment(), 
respectively.

For information about the Gblocks program, see ????

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:


=head1 AUTHOR - Robson Francisco de Souza

Email rfsouza-at-cecm.usp.br

=head1 CONTRIBUTORS

This code was heavyly based on Jason Stajich's
Bio::Tools::Run::Phylo::PAML::Yn00.pm module. 

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Tools::Run::Alignment::Gblocks;

use IO::String;
use Bio::Root::Root;
use Bio::AlignIO;
use Bio::Tools::Run::WrapperBase;
use strict;

our @ISA = qw(Bio::Root::Root Bio::Tools::Run::WrapperBase);

# Setting default values and allowed parameters
our $PROGRAMNAME = 'Gblocks';
our $PROGRAMDIR  = '/usr/local/bin';

our %VALID_PARAMETERS = 
    (
     'min_init_block_length'       => 'b0', 'b0' => 'b0',
     'min_number_seqs_conserved'   => 'b1', 'b1' => 'b1',
     'min_number_seqs_flanking'    => 'b2', 'b2' => 'b2',
     'max_contiguous_nonconserved' => 'b3', 'b3' => 'b3',
     'min_block_length'            => 'b4', 'b4' => 'b4',
     'allowed_gaps'                => 'b5', 'b5' => 'b5',
     'similarity_matrix'           => 'b6', 'b6' => 'b6',
     'file_extension'              => 'e',  'e'  => 'e',
     'save_masked_alignment'       => 'k',  'k'  => 'k',
     'save_nonconserved_blocks'    => 'n',  'n'  => 'n',
     'save_results_and_parameters' => 'p',  'p'  => 'p',
     'save_selected_blocks'        => 's',  's'  => 's',
     'type'                        => 't',  't'  => 't',
     'save_ungapped_alignment'     => 'u',  'u'  => 'u',
     'characters_per_line'         => 'v',  'v'  => 'v',
     );

our %DEFAULT_PARAMETERS = 
    (
     'allowed_gaps'                => undef,
     'characters_per_line'         => undef,
     'file_extension'              => undef,
     'max_contiguous_nonconserved' => undef,
     'min_block_length'            => undef,
     'min_init_block_length'       => undef,
     'min_number_seqs_conserved'   => undef,
     'min_number_seqs_flanking'    => undef,
     'save_selected_blocks'        => 'y',
     'save_results_and_parameters' => 'Short Text',
     'similarity_matrix'           => undef,
     'type'                        => undef,
     'save_masked_alignment'       => undef,
     'save_nonconserved_blocks'    => undef,
     'save_ungapped_alignment'     => undef,
     );

=head1 Object constructor

=head2 new

 Title   : new
 Usage   : my $obj = new Bio::Tools::Run::Align::Gblocks();
 Function: Builds a new Bio::Tools::Run::Align::Gblocks object
 Returns : Bio::Tools::Run::Align::Gblocks
 Args    : named arguments. Available options are
           # General options
           -verbose          => supported by most Bioperl modules, when set
                                to true (1) this flag activates printing of
                                standard output messages from GBlocks
           -executable       => set full path to Gblocks executable
           -save_tempfiles   => boolean to save the generated tempfiles and
                                NOT cleanup after onesself (default FALSE)

           # Specifying an alignent:
           -alignment        => loaded alignment (Bio::Align::AlignI object)

           # or
           -alignment_file   => load alignment from this file
           -alignment_format => input file format (Bio::AlignIO compatible)

 Note    : option -alignment_format has no effect without -alignment_file

=cut

sub new {
  my($class,@args) = @_;

  my $self = $class->SUPER::new(@args);
  $self->set_default_parameters();
  $self->program_name($PROGRAMNAME);
  $self->program_dir($PROGRAMDIR);
  $self->{'_arguments'} = undef;
  $self->{'_tempdir'} = 0;

  my ($st,$prog,$quiet,$aln,$alnfile,$format,$tempdir) = 
      $self->_rearrange([qw(SAVE_TEMPFILES EXECUTABLE QUIET ALIGNMENT
			    ALIGNMENT_FILE ALIGNMENT_FORMAT TEMPDIR
			    )], @args);

  $format = 'phylip' unless (defined $format);
  $self->{'_tempdir'} = $tempdir if (defined $tempdir);
  defined $prog && $self->executable($prog);
  defined $quiet && $self->quiet($quiet);
  defined $st  && $self->save_tempfiles($st);
  $self->alignment($aln) if (defined $aln);
  $self->alignment_file($alnfile,$format) 
      if (defined $alnfile && !defined $aln);

  return $self;
}

=head1 Class specific methods

=head2 alignment

 Title   : alignment
 Usage   : $self->alignment($name);
 Function: Get/Set the input alignment object
 Returns : L<Bio::Align::AlignI>
 Args    : a Bio::Align::AlignI object

 Comment : see L<alignment_file> for info on how to
           use an alignment file, instead of an object,

=cut

sub alignment {
   my ($self,$obj) = @_;

   if (defined $obj) {
       if (ref($obj) && $obj->isa("Bio::Align::AlignI")) {
	   $self->{'_alignment'} = $obj;
       } else {
	   $self->throw(qq{Use method alignment_file to load alignment from disk}) if (-f $obj);
	   $self->throw(qq{Expecting Bio::Align::AlignI object, instead of $obj});
       }
   }

   return $self->{'_alignment'};
}

=head2 alignment_file

 Title   : alignment_file
 Usage   : $self->alignment_file($name, $format);
 Function: Get/Set the input alignment file
 Returns : file name or undef
 Args    : file name and (optional) Bio::AlignIO 
           format, i.e. driver name.

 Comment : if the file is not Phylip formatted, it
           will be converted to Phylip transparently.
           See L<alignment> for info on how to use
           an alignment object, instead of a file.

=cut

sub alignment_file {
   my ($self,$file,$format) = @_;

   if (defined $file && -e $file) {
       if (defined $format && $format !~ /^fasta$/i) {
	   $self->warn("Loading non-PHYLIP formatted input alignment");
	   my $in = Bio::AlignIO->new(-file=>$file, -format=>$format);
	   $self->{'_alignment'} = $in->next_aln;
	   $in->close(); undef($in);
       }
       $self->{'_alignment_file'} = $file;
   }

   return $self->{'_alignment_file'};
}

=head2 get_parameter

 Title   : get_parameter
 Usage   : my %params = $self->get_parameter($name);
 Function: Get the value set for a parameter
 Returns : current parameter value
 Args    : parameter name or undef (if not set)

=cut

sub get_parameter {
   my ($self,$name) = @_;
   $self->warn(qq{Requesting value for unknown parameter $name})
       unless (exists $VALID_PARAMETERS{$name});
   $name = $VALID_PARAMETERS{$name};
   return undef unless (exists $self->{'_parameters'}{$name});
   return $self->{'_parameters'}{$name};
}

=head2 get_all_parameters

 Title   : get_all_parameters
 Usage   : my %params = $self->get_all_parameters();
 Function: returns the list of parameters as a hash
 Returns : associative array keyed on parameter names
 Args    : none

=cut

sub get_all_parameters {
   my ($self) = @_;
   # we're returning a copy of this
   return %{ $self->{'_parameters'} };
}

=head2 set_parameters

 Title   : set_parameters
 Usage   : $wrapper->set_parameters(%hash);
 Function: Set Gblocks parameters. 

           This module does not implements all Gblocks options
           and supported parameters are validated using the
           class variable %VALID_PARAMETERS. To pass arbitrary
           parameters to Gblocks, one may turn off parameter
           checking using

             $wrapper->no_param_checks(1).

 Returns : boolean if all parameters are valid
 Args    : Hash of parameters (keys) and new values

 See also: 
           L<no_param_checks()>

           Gblocks documentation for command line options 
           and the contents of %VALID_PARAMETERS, e.g. using this code:

           use Bio::Tools::Run::Align::Gblocks;
           my %a = %Bio::Tools::Run::Align::Gblocks::VALID_PARAMETERS;
           print join("\n",keys %a),"\n";

=cut

sub set_parameters {
   my ($self,%hash) = @_;

   my $ret = 1;
   while (my ($param,$value) = each %hash) {
       unless ($self->no_param_checks || exists $VALID_PARAMETERS{$param}) {
	   $self->warn(qq{Gblocks: Unknown parameter $param will not set unless you force by setting no_param_checks to true});
	   $ret = 0;
	   next;
       }
       # Turning parameters into command line switches
       $param = $VALID_PARAMETERS{$param};
       $self->{'_parameters'}->{$param} = $value;
   }

   $self->{'_arguments'} = undef;
   return $ret;
}

=head2 set_default_parameters

 Title   : set_default_parameters
 Usage   : $wrapper->set_default_parameters()
 Function: (Re)set all parameters to default values
 Returns : nothing
 Args    : boolean: keep existing parameter values

=cut

sub set_default_parameters {
    my ($self,$keepold) = @_;
    $keepold = 0 unless defined $keepold;

    while( my ($param,$val) = each %DEFAULT_PARAMETERS ) {
	# Turning parameters into command line switches
	$param = $VALID_PARAMETERS{$param};
	# skip if we want to keep old values and it is already set
	next if(defined $self->{'_parameters'}->{$param} && $keepold);
	$self->{'_parameters'}->{$param} = $val;
    }
    $self->{'_arguments'} = undef;
}

=head1 Bio::Tools::Run::WrapperBase methods implemented by this class

=head2 arguments

 Title   : arguments
 Usage   : $obj->arguments($newval)
 Function: Commandline parameters
 Returns : value of arguments
 Args    : string (optional)

 Comments: Although setting Gblocks commad line arguments
           by assigning a string to $newval is left as
           an option to the user, the Gblocks command line
           is usually built from parameters set by using
            L<set_parameters> or L<set_default_parameters>.

=cut

sub arguments {
    my ($self,$args) =  @_;

    unless (defined $args || defined $self->{'_arguments'}) {
	$args = '';
	my %param = $self->get_all_parameters;
	foreach my $arg (keys %param) {
	    next unless (defined $param{$arg});
	    $args .= " -${arg}=".$param{$arg};
	}
    }

    $self->{'_arguments'} = $args;
    return $self->{'_arguments'};
}

=head2 run

 Title   : run
 Usage   : ($rc,$data) = $wrapper->run($aln);
 Function: Actually running Gblocks
 Returns : run status (1 if ok, 0 if error) and
           hash reference with keys
            selected_alignment => reduced alignment with
                                  selected blocks
            ungapped_alignment => reduced alignment with
                                  selected blocks
            results => reference to a hash describing the
                       results of the run.

 Args    : (optional) L<Bio::Align::AlignI>
           object or alignment/alignment path filename

 Comments: the argument makes the method ignore
           any alignment object or file set
           using methods alignment() or
           alignment_file().
           
           The user can also use as argument the path
           to a file listing all alignments to analyze,
           which is an option avaliable in Gblocks.

=cut

sub run {
    my ($self,$aln) = @_;

    # Preparing alignment file
    $aln = $self->alignment unless (defined $aln);
    my $alnfile = undef;
    if (!defined $aln) {
	$alnfile = $self->alignment_file;
    } elsif (!ref($aln) && -e $aln) {
	$alnfile = $aln;
    } elsif ($aln eq "-") {
	my $alignFH;
	($alignFH,$alnfile) = $self->io->tempfile(UNLINK => ($self->save_tempfiles ? 0 : 1), DIR => '.');
	while (<>) { print $alignFH $_ };
	close($alignFH); undef($alignFH);
    } elsif (ref($aln) && $aln->isa('Bio::Align::AlignI')) {
	my $alignFH;
	($alignFH,$alnfile) =
	    $self->io->tempfile(UNLINK => ($self->save_tempfiles ? 0 : 1),
				DIR    => '.');
	my $alnout = Bio::AlignIO->new(-format => 'fasta',
				       -fh     => $alignFH,
				       -displayname_flat => 1);
	$alnout->write_aln($aln);
	$alnout->close(); undef($alnout);
	close($alignFH); undef($alignFH);
    } else {
	$self->throw(qq{Unknown alignment argument $aln});
    }
    $self->throw(qq{No valid alignment loaded when trying to run Gblocks!})
	unless (defined $alnfile && -e $alnfile);

    # Preparing command line
    my $gblocks = $self->executable();
    $self->throw('Unable to find executable '.
		 $self->program_name.' for Gblocks') unless (defined $gblocks);
    my $cmd = "$gblocks $alnfile ".$self->arguments;
    $cmd .= ' > /dev/null' unless ($self->verbose);

    # Run!
    my ($rc,$data) = (0,{});
    eval { $rc = system($cmd) };

    # Processing Gblocks output
    my $ext      = $self->get_parameter('e') || '-gb';
    my $save_aln = $self->get_parameter('s');
    eval {
	# Loading selected alignment columns
	unless (defined $save_aln && $save_aln eq 'n') {
	    open(FASTA,"<${alnfile}${ext}") || $self->throw("Could not open Gblocks Fasta output file ${alnfile}${ext}");
	    my $fasta = "";
	    while (<FASTA>) {
		s/ //g unless /^>/;
		$fasta .= $_;
	    }
	    close(FASTA);
	    my $str = IO::String->new($fasta);
	    my $aio = Bio::AlignIO->new(-fh=>$str, -format=>'fasta');
	    $data->{'selected_alignment'} = $aio->next_aln;
	    $aio->close; undef($aio);
	}
	# Loading ungapped alignment
	my $save_ungapped = $self->get_parameter('u');
	if (defined $save_ungapped && $save_ungapped eq 'y') {
	    my $aio = Bio::AlignIO->new(-file=>"${alnfile}${ext}--",
					-format=>'fasta');
	    $data->{'ungapped_alignment'} = $aio->next_aln;
	}
	# Loading run parameters and results
	$data->{'results'} = $self->_load_gblocks_results($alnfile);
    };
    $self->throw("Error loading Gblocks output for $alnfile:\n$@") if ($@);

    # Cleaning up: removing Gblocks output files 
    # and temporary files and/or directories
    if ($self->save_tempfiles) {
	$data->{'results'}{'selected_alignment_file'} = "${alnfile}${ext}"
	    if ( -e "${alnfile}${ext}");
    } else {
	unlink("${alnfile}${ext}")      if ( -e "${alnfile}${ext}");
	unlink("${alnfile}${ext}--")    if ( -e "${alnfile}${ext}--");
	unlink("${alnfile}${ext}.txts") if ( -e "${alnfile}${ext}.txts");
	unlink("${alnfile}${ext}.txt")  if ( -e "${alnfile}${ext}.txt");
	unlink("${alnfile}${ext}.htm")  if ( -e "${alnfile}${ext}.htm");
	unlink("${alnfile}${ext}Mask")  if ( -e "${alnfile}${ext}Mask");
	unlink("${alnfile}${ext}PS")    if ( -e "${alnfile}${ext}PS");
	unlink("${alnfile}${ext}Comp")  if ( -e "${alnfile}${ext}Comp");
	unlink("${alnfile}.seq")        if ( -e "${alnfile}.seq");
	$self->cleanup();
    }

    return ($rc,$data);
}

=head2 program_name

 Title   : program_name
 Usage   : $wrapper->program_name($newvalue)
 Function: Get/Set the program name
 Returns : string (default: Gblocks)
 Args    : string

=cut

sub program_name {
    my ($self) = shift;
    $self->{'_program_name'} = shift if (scalar @_);
    return $self->{'_program_name'};
}

=head2 program_dir

 Title   : program_dir
 Usage   : $wrapper->program_dir()
 Function: Get/Set the directory where the program is
 Returns : string (default: /usr/local/bin)
 Args    : string 

=cut

sub program_dir {
    my ($self) = shift;
    $self->{'_program_dir'} = shift if (scalar @_);
    return $self->{'_program_dir'};
}

=head1 Bio::Tools::Run::WrapperBase inherited methods

=head2 error_string

 Title   : error_string
 Usage   : $obj->error_string($newval)
 Function: Where the output from the last analysis run is stored.
 Returns : value of error_string
 Args    : newvalue (optional)

=head2 no_param_checks

 Title   : no_param_checks
 Usage   : $obj->no_param_checks($newval)
 Function: Boolean flag as to whether or not we should
           trust the sanity checks for parameter values
 Returns : value of no_param_checks
 Args    : newvalue (optional)

=head2 save_tempfiles

 Title   : save_tempfiles
 Usage   : $obj->save_tempfiles($newval)
 Function:
 Returns : value of save_tempfiles
 Args    : newvalue (optional)

=head2 outfile_name

 Title   : outfile_name
 Usage   : my $outfile = $wrapper->outfile_name();
 Function: Get/Set the name of the output file for this run
           (if you wanted to do something special)
 Returns : string
 Args    : [optional] string to set value to

=head2 tempdir

 Title   : tempdir
 Usage   : my $tmpdir = $self->tempdir();
 Function: Retrieve a temporary directory name (which is created)
 Returns : string which is the name of the temporary directory
 Args    : none

=head2 cleanup

 Title   : cleanup
 Usage   : $wrapper->cleanup();
 Function: Will cleanup the tempdir directory
 Returns : none
 Args    : none

=head2 io

 Title   : io
 Usage   : $obj->io($newval)
 Function: Gets a Bio::Root::IO object
 Returns : Bio::Root::IO object
 Args    : none

=head2 version

 Title   : version
 Usage   : $version = $wrapper->version()
 Function: Returns the program version (if available)
 Returns : string representing version of the program
 Args    : [Optional] value to (re)set version string

=head2 executable

 Title   : executable
 Usage   : my $exe = $factory->executable();
 Function: Finds the full path to the executable
 Returns : string representing the full path to the exe
 Args    : [optional] name of executable to set path to
           [optional] boolean flag whether or not warn when exe is not found

=head2 program_path

 Title   : program_path
 Usage   : my $path = $factory->program_path();
 Function: Builds path for executable
 Returns : string representing the full path to the exe
 Args    : none

=head1 Internal methods

=head2 _load_gblocks_stats

 Title   : _load_gblocks_stats
 Usage   : $wrapper->_load_gblocks_stats($name);
 Function: Load Gblocks detailed block selection
           information from one of these files: 

           ${name}${ext}.txts
           ${name}${ext}.txt
           ${name}${ext}.html 

 Returns : hash reference
 Args    : input alignment file name

=cut

sub _load_gblocks_results {
    my ($self,$aln) = @_;
    my $ext = $self->get_parameter("e");
    $ext = "-gb" unless (defined $ext);

    my $resfile = "${aln}${ext}.txts";
    $resfile = "${aln}${ext}.txt" if (! -e $resfile);
    $resfile = "${aln}${ext}.htm" if (! -e $resfile);

    open(STAT,"<$resfile") || $self->throw("Could not open result file $resfile for alignment $aln");
    my $stat = {};
    while (<STAT>) {
	chomp;
	next if (/^\s*$/);

	/Alignment assumed to be:\s+(\<b\>)*(\w+)(\<b\>)*/ && do {
	    $stat->{'type'} = $2;
	};

	/New number of positions: (\d+)/ && do {
	    $stat->{'new_length'} = $1;
	};

	/Minimum Number Of Sequences For A Conserved Position: (\d+)/ && do {
	    $stat->{'min_number_seqs_conserved'} = $1;
	};

	/Minimum Number Of Sequences For A Flanking Position: (\d+)/ && do {
	    $stat->{'min_number_seqs_flanking'} = $1;
	};

	/Maximum Number Of Contiguous Nonconserved Positions: (\d+)/ && do {
	    $stat->{'max_contiguous_nonconserved'} = $1;
	};

	/Minimum Length Of An Initial Block: (\d+)/ && do {
	    $stat->{'min_init_block_length'} = $1;
	};

	/Minimum Length Of A Block: (\d+)/ && do {
	    $stat->{'min_block_length'} = $1;
	};

	/Allowed Gap Positions: (.+)/ && do {
	    $stat->{'allowed_gaps'} = $1;
	};

	/Use Similarity Matrices: (\d+)/ && do {
	    $stat->{'similarity_matrix'} = $1;
	};

	/Flanks:\s+\[(.+)\]/ && do {
	    my $flanks = $1;
	    my @flanks = map { [ split(/\s+/,$_) ] } split(/\]\s+\[/,$flanks);
	    $stat->{'flanks'} = \@flanks;
	};


	/\(\d+\% of the original (\d+) positions\)/ && do {
	    $stat->{'original_length'} = $1;
	};
    }
    close(STAT);

    return $stat;
}

1;
