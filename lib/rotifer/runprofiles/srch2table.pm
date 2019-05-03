# Default runprofiles postprocessor
sub process {
    my ($config, $target, $in, $out, $err) = @_;

    # Remove empty input files
    if ( -z $in ) {
	open (ERR,">>$err") || die "Could not open standard error file $err";
	print ERR "File $in was empty (no results): removed!";
	close(ERR);
	unlink($in);
	return -1;
    } 

    # Get target name
    my $name = $in;
    if (defined $target && ref $target) {
	if ($target->can("domain")) {
	    $name = $target->domain;
	} elsif ($target->can("name")) {
	    $name = $target->name;
	} elsif ($target->can("id")) {
	    $name = $target->id;
	}
    }

    # Build commands
    my @command = ();
  SWITCH: {
      $_ = $in; # Simplify regexps below

      # BLAST
      /\.(psi|t)?blast(all|pgp|n|p|x)?$/ && do {
	  push(@command, "blast2table -r $name -t ".$config->threshold." -e ".$config->pcut." $in > $out 2> $err");
	  last SWITCH;
      };

      # RPS-BLAST
      /\.rpsblast$/ && do {
	  my $clean = /(Cdd|Pfam|Smart)/i ? 'cdd' : 'profiledb';
	  push(@command, "blast2table -s -c $clean -t ".$config->threshold." -e ".$config->pcut." $in > $out 2> $err");
	  last SWITCH;
      };

      # HMMER
      /\.(hmmsearch|jackhmmer)$/ && do {
	  push(@command, "hmmer2table -r $name -c -e ".$config->pcut." $in > $out 2> $err");
	  last SWITCH;
      };
      /\.hmmscan$/ && do {
	  push(@command, "hmmer2table -c -s -e ".$config->pcut." $in > $out 2> $err");
	  last SWITCH;
      };

      open (ERR,">>$err") || die "Could not open standard error file $err";
      print ERR "$in: unknown file format";
      close(ERR);
      exit 1;
    };
    push(@command, "gzip $in 2>> $err");
    unlink("${in}.gz") if ( -f "${in}.gz" );

    # Run
    my $stat = 0;
    foreach my $cmd (@command) {
	($stat) = (system("$cmd") >> 8); # 0 if ok
	if ($stat) { # If system() call failed...
	    print STDERR "Error running postprocessor srch2table for batch $out (shell command: $cmd)";
	    exit 2;
	}
    }

    # remove empty tables
    unlink($out) if ( -z $out );
    unlink($err) if ( -z $err );

    return $stat;
}

1; # To make perl happy
