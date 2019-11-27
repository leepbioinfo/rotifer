# Default runprofiles postprocessor
sub _find_executable {
    use File::Which qw(which);
    my $executable = undef;
    foreach my $name (@_) {
        next unless (defined $name && length $name);
        foreach my $exe (split(/\s*\,\s*/,$name), which($name)) {
            if ( -x $exe ) {
                $executable = $exe;
                last;
            }
        }
    }
    return $executable;
}

####
# Check for overlaps in BOTH the query and hit sequences
# This subroutine supports tolerance of overlaps based on percentage
# (float, e.g. 1.0) or number of residues (integer, e.g. 1).
sub same_hsp {
    my ($config, $array1, $array2) = @_;
    # saccver: 0, qaccver: 1, sstart: 2, send: 3, evalue: 4, qcovhsp: 5, qstart: 6, qend: 7, bitscore: 8, slen: 9

    my $overlaps = 0; # 0 Means $array1 and $array2 are not the same
    my $maxoverlap = $config->postprocessor_param->{'maxoverlap'} || 0.4;
    foreach my $i (2,6) {
        my $max_start = $array1->[$i]   > $array2->[$i]   ? $array1->[$i]   : $array2->[$i];
	my $min_end   = $array1->[$i+1] < $array2->[$i+1] ? $array1->[$i+1] : $array2->[$i+1];
        my $overlap_length = $min_end - $max_start + 1;
	#print join(" ","Check $i:",$maxoverlap,$overlap_length,@$array1,@$array2),"\n";
	if ($maxoverlap =~ /\d*\.\d+/ || $maxoverlap =~ /\d+e-\d+/) { # Compare both HSPs when using percentage!
            my $max_overlap1 = $maxoverlap * ($array1->[$i+1] - $array1->[$i] + 1);
            my $max_overlap2 = $maxoverlap * ($array2->[$i+1] - $array2->[$i] + 1);
            $overlaps++ if ($overlap_length > $max_overlap1 && $overlap_length > $max_overlap2);
        } else { # Cutoff set in number of residues: compare length of overlapping region to cutoff
            $overlaps++ if ($overlap_length > $maxoverlap);
        }
    }

    return $overlaps == 2 ? 1 : 0; # 0 => != HSPs, 1 => hit overlap, 2 => query and hit overlap
}

sub process {
    my ($config, $target, $in, $out, $err) = @_;
    my $blast2table = _find_executable($config->program_path->{"blast2table"}) || which("blast2table") || "blast2table";
    my $hmmer2table = _find_executable($config->program_path->{"hmmer2table"}) || which("hmmer2table") || "hmmer2table";

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
	  #push(@command, "$blast2table -r $name -t ".$config->threshold." -e ".$config->pcut." $in > $out 2> $err");
	  open(IN,"<$in")   || die "Could not open search results $in";
	  my $iteration = 1;
	  my %index  = ();
	  my @parsed = ();
	  while (<IN>) {
	      # saccver: 0, qaccver: 1, sstart: 2, send: 3, evalue: 4, qcovhsp: 5, qstart: 6, qend: 7, bitscore: 8, slen: 9
	      chomp;
	      $H{i} = $1 if (/# Iteration:\s+(\d+)/);
	      next if (/^\#|^\s*$|^Search has CONVERGED/);
	      my @F = split("\t");
	      @F[8..11] = ($H{i},@F[8..9],$F[1]);
	      $F[1] =~ s/(\_\S+)?\.\d+$//;
	      $F[5] = sprintf("%.3f",$F[5]/100);

	      # Check for redundant HSPs, i.e. HSPs that overlap both in the query and the hit
	      my $id = "$F[0]:$F[1]";
	      my $is_new = 1;
	      if (exists($H{i}) && exists($index{$id})) { # Second HSP of a query/subject pair from PSIBLAST?
		  foreach my $i (@{ $index{$id} }) {
		      #print join(" ","Check",$H{i},$i),"\n";
		      if (same_hsp($config, $parsed[$i], \@F)) {
			  if ($parsed[$i]->[4] < $config->threshold) {
			      # If we have already found the first HSP lower than the inclusion threshold, we
			      # will select and store the widest coordinates but keep the E-value and iteration
			      foreach my $j (2,6) { # 2: hit start, 6: query start
				  my $old = $parsed[$i]->[$j+1] - $parsed[$i]->[$j];
				  my $new = $F[j] - $F[j];
				  if ($new > $old) {
				      $parsed[$i]->[$j+1] = $F[${j+1}];
				      $parsed[$i]->[$j]   = $F[$j];
				  }
			      }
			  } elsif ($F[4] < $parsed[$i]->[4]) { # 4: evalue
			      # Replace a previously loaded HSP with one with lower E-value
			      $parsed[$i] = [ @F ]; # Copy
			  }
			  $is_new = 0; # If we made it to this line, then we already loaded this HSP
			  last;        # Abort inspecting loaded HSPs when the first equivalent HSP is found
		      } # if (same_hsp($config, $parsed[$i], $data))
		  } # foreach my $i (@{ $index{$id} })
	      } # if (exists($H{i}) && exists($index{$id}))

	      # Add new HSP
	      if ($is_new) {
		  my $copy = [ @F ];
		  #&swap($copy) if ($config->swap);
		  push(@parsed, $copy); # Copy each HSP's data to this array
		  push(@{ $index{$id} }, $#parsed);
	      }
	  }
	  close(IN);

	  # Print!
	  if (scalar @parsed) {
	      open(OUT,">$out") || die "Could not open processed output file $out for printing";
	      print OUT join("\n",map { join("\t",@$_) } @parsed),"\n";
	      close(OUT);
	  }
	  last SWITCH;
      };

      # RPS-BLAST
      /\.rpsblast$/ && do {
	  my $clean = /(Cdd|Pfam|Smart)/i ? 'cdd' : 'profiledb';
	  push(@command, "$blast2table -s -c $clean -t ".$config->threshold." -e ".$config->pcut." $in > $out 2> $err");
	  last SWITCH;
      };

      # HMMER
      /\.(hmmsearch|jackhmmer)$/ && do {
	  push(@command, "$hmmer2table -r $name -c -e ".$config->pcut." $in > $out 2> $err");
	  last SWITCH;
      };
      /\.hmmscan$/ && do {
	  push(@command, "$hmmer2table -c -s -e ".$config->pcut." $in > $out 2> $err");
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
