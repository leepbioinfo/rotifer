# Default runprofiles postprocessor
#
sub process {
    my ($config, $target, $in, $out, $err) = @_;

    push(@command, "gzip $in 2>> $err");

    # Run
    my $stat = 0;
    foreach my $cmd (@command) {
        ($stat) = (system("$cmd") >> 8); # 0 if ok
        if ($stat) { # Cannot use system BLAST database
            print STDERR "Error running postprocessor srch2table for batch $out (shell command: $cmd)";
            exit 2;
        }
    }

    return $stat;
}

1; # To make perl happy
