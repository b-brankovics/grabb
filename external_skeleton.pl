#!/usr/bin/perl -w
use strict;

# A perl script that allows you to use an assembler that is not specified in the GRAbB source code

# Get all the parameters
my ($reads, $paired, $parameters1, $parameters2, $outfile, $format) = @ARGV;
# Create an array containing the read files
my @reads = split / /, $reads;
my @param1 = split / /, $parameters1;
my @param2 = split / /, $parameters2;

# If there were no extra arguments passed for the assembler at the invocation of GRAbB
#   then the value passed to this script is "-" ($param1[0] and/or $param2[0])

# For mapping-assemblers a reference is needed
#my $ref = &get_reference();  # For reference based assembly uncomment this line


# The program to be used for assembly
my $assmebler = ""; # ADD

# Run the assembly
`$assmebler @reads`; # ADD

# Specifiy the expected result file
my $result = ""; # ADD

# Copy the result file to the expected position
`perl -ne 'if (/^>/) {print} else {tr/acgtwmrskybvdh/ACGTWMRSKYBVDH/; print}' $result >$outfile`;

sub get_reference {
    # This subroutine finds the reference file for the current thread
    my $ref;
    my $pwd = $ENV{"PWD"};
    $pwd =~ /(thread_\d+)$/;
    my $thread = $1;
    $ref = "../../$thread/bait.fas";
    return $ref;
}
