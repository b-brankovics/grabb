#!/usr/bin/perl -w
use strict;

############################################################
# A perl script that allows you to use SPades assembler    #
# that is not specified in the GRAbB executable            #
############################################################

#                      +----+
#                      |    |
#                      |    |
#                      |    |
#                      |    |
#                    \-+    +-/
#                     \      /
#                      \    /
#                       \  /
#                        \/

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Specify the executable for SPAdes
my $assmebler = "spades.py"; # ADD executable here!

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#                        /\
#                       /  \
#                      /    \
#                     /      \
#                    /-+    +-\
#                      |    |
#                      |    |
#                      |    |
#                      |    |
#                      +----+



unless ($assmebler) {
    print "No assembler command was specified!\n";
    die "No assembler command was specified! Please, update the external skeleton file!\n";
}


if (scalar(@ARGV) == 0 || grep{/^-?-h(elp)?$/} @ARGV) {
    print "USAGE:\n\tperl $0 " .
	'"<readfile1> <readfile2> ..." ' .
	"<paired|single> " . 
	'"<assembly paramters>" ' .
	'"<input parameters>" ' .
	"<assembly.fas>\n";

    print "\nDescription:\n",
    "\tA perl script that allows you to use SPades assembler\n",
    "\tthat is not specified in the GRAbB executable\n\n";

    print "TODO:\n",
    "\tAdd the desired 'spades.py' file to path ('spades.py' should an exutable command)\n",
    "\t\tif 'spades.py' is in the path, then check steps 3. and ! below\n",
    "\t\t3. Run GRAbB.pl by adding the following:\n\t\t\t--assembler external " . $ENV{"PWD"} ."/$0\n",
    "\t\t!  Specify assembly parameters using '--arg1' option of GRAbB.pl. For example: --arg1 -k 31,61,91\n",
    "\tOR do the following steps:\n",
    "\t1. Open this file in a text editor and specify the excutable of SPAdes\n",
    "\t\tIf the SPAdes command is available as 'spades-3.10.0.py', then\n",
    "\t\tModify line 22:\n\t\t" . 'my $assmebler = "spades.py"; # ADD executable here!' . "\n",
    "\t\t                  ||\n",
    "\t\t                  \\/\n",
    "\t\t" . 'my $assmebler = "spades-3.10.0.py"; # ADD executable here!' . "\n\n",
    "\t\tOR if the SPAdes executable's path is '~/bin/spades.py', then\n",
    "\t\tModify line 22:\n\t\t" . 'my $assmebler = "spades.py"; # ADD executable here!' . "\n",
    "\t\t                  ||\n",
    "\t\t                  \\/\n",
    "\t\t" . 'my $assmebler = "~/bin/spades.py"; # ADD executable here!' . "\n",
    "\n",
    "\t2. Save the modified file\n",
    "\t3. Run GRAbB.pl by adding the following:\n\t\t--assembler external " . $ENV{"PWD"} ."/$0\n",
    "\t!  Specify assembly parameters using '--arg1' option of GRAbB.pl. For example: --arg1 -k 31,61,91\n",
    "\n\n";
    
    print "Example parameters (This is only for testing purposes):\n",
    "\t" . '"' . "<readfile1> <readfile2> ..." . '"'. "\n\t\t" . '"'. "readpool1.fastq readpool2.fastq" . '"'. "\n",
    "\t<paired|single>\n\t\tpaired\n",
    "\t" . '"<assembly parameters>"' . "\n\t\t" . '"' . "-k 31,61,91" . '"' . "\n",
    "\t" . '"<input parameters>" '   . "\n\t\t" . '"' . "--pe1-1 readpool1.fastq --pe1-2 readpool2.fastq" . '" ' . "\n",
    "\t" . "<assembly.fas>"          . "\n\t\t" . "assembly.fas\n";
    exit;
}

# `perl ../../assembler.pl "@$readpool_ref" $arg_single "$extra1" "$extra2" $$out_file_ref $$format_ref`;   # Assemble
# Get all the parameters
my ($reads, $paired, $parameters1, $parameters2, $outfile, $format) = @ARGV;
# Create an array containing the read files
my (@reads, @param1, @param2);

# If there were no extra arguments passed for the assembler at the invocation of GRAbB
#   then the value passed to this script is "-" ($param1[0] and/or $param2[0])
unless ($parameters1 eq "-") {
    @param1 = split / /, $parameters1;
}
unless ($parameters2 eq "-") {
    @param2 = split / /, $parameters2;
    $paired = "skip";
}

@reads = split / /, $reads;

# Build 'read_line' for single or paired options
my $read_line = "";
if ($paired eq "paired") {
    if (scalar(@reads) % 2 == 0) {
	my $lib = 0;
	while ($lib < scalar(@reads) / 2) {
	    my $f = $lib * 2;
	    my $r = $f + 1;
	    $lib++;
	    $read_line .= " --pe" . $lib . "-1 $reads[$f] --pe" . $lib . "-2 $reads[$r]";
	}
    } else {
	print "Uneven number of read files, while using 'paired mode'!\n";
	die "Uneven number of read files, while using 'paired mode'!\n";
    }
} elsif ($paired eq "single") {
    my $lib;
    for (@reads) {
	$lib++;
	$read_line .= " --s$lib $_";
    }
} else {
    # Overwrite default read input info
    $read_line = join(" ", @param2);
}



# For mapping-assemblers a reference is needed
my $ref = &get_reference();  # For reference based assembly uncomment this line


# Specifiy the expected result file
my $folder = "spades";
if (grep{$_ eq "-o"} @param1) {
    # If different output folder is specified then change result folder
    # Get the index of the folder name
    my $i;
    for (@param1) {
	$i++;
	if ($_ eq "-o") {
	    last;
	}
    }
    $folder = $param1[$i];
} else {
    push @param1, "-o spades";
}
my $file = "contigs.fasta";
my $result = $folder . "/" . $file;

unless ($result) {
    print "No result file was specified!\n";
    die "No result file was specified! Please, update the external skeleton file!\n";
}

# Default assembly parameters if none is specified
my $assembly_param = "-k auto"; # my typical setting for 100bp-long reads: '-k 31,61,91'
my $param_defined = grep{/^-k/} @param1;
if (scalar(@param1) == 0 || ! $param_defined) {
    push @param1, $assembly_param;
}


# Run the assembly
`$assmebler $read_line @param1`;
#print "$assmebler $read_line @param1\n";

my $gfa = "spades/assembly_graph_with_scaffolds.gfa";

my $select = "selected.gfa";

`gfa_steps.pl $ref $gfa $select $outfile`;

# Copy the result file to the expected position
#`perl -ne 'if (/^>/) {print} else {tr/acgtwmrskybvdh/ACGTWMRSKYBVDH/; print}' $result >$outfile`;

sub get_reference {
    # This subroutine finds the reference file for the current thread
    my $ref;
    my $pwd = $ENV{"PWD"};
    $pwd =~ /(thread_\d+)$/;
    my $thread = $1;
    $ref = "../../$thread/reference.fas";
    return $ref;
}
