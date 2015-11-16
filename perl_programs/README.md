# Perl programs

These are programs that use only perl and basic unix commands and do
not call any programs that need to be separately installed.

- kmer\_bait\.pl

    + It is an alternative to using mirabait for exact 31-mer
      matching.
	+ It is less efficient than mirabait. Therefore it is recommended
      to install mirabait (not v 4.0.2).

- create\_readpool\.pl

    + It is an alternative to using seqtk to extract reads that
      matched to the bait
	+ It is less efficient than seqtk. Therefore it is recommended
      to install seqtk.

## Add these programs to the GRAbB\.pl source code

- This can be done using the **configure_GRAbB.pl** script (this adds
  either these programs or the recommended programs for these steps if
  they are found by the script and are working as expected) or manually:

1. Make the files executable

        chmod +x <file_name>

2. Add the executable to the path or to your __bin__ directory

        cp <file_name> /home/user/bin/.

3. Adjust the source code of GRAbB\.pl

    - kmer\_bait\.pl

        Replace the following line:

            my $bait_cmd =    "mirabait";       # The command to invoke the baiting program                    #	!!!!

        with:

            my $bait_cmd =    "kmer_bait.pl";       # The command to invoke the baiting program                    #	!!!!

    - create\_readpool\.pl

        Replace the following line:

            my $collect_cmd = "seqtk subseq";   # The command to invoke the read collecting program            #	!!!!

        with:

            my $collect_cmd = "create_readpool.pl";   # The command to invoke the read collecting program            #	!!!!
