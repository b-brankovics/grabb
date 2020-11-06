#!/usr/bin/perl -w
use strict;
use autodie;

# This is a perl alternative for using mirabait to bait reads based on protein sequences using kraken2
# Invocation by GRAbB.pl:
#    $bait_cmd -t txt $$bait_ref $read $mira_temp\_$i 2>&1 >>mirabait.log

# Input:
#        "-t"
#        "txt"
#        Reference (Fasta)
#        Read file (Fasta or Fastq)
#        prefix (prefix for the output file: prefix.txt)

# kraken executables
my $kbcmd = "kraken2-build";
my $kcmd = "kraken2";

# Test correct invocation
my $usage = "Usage\n\t$0 -t txt <reference.fas> <read.fastq|read.fasta> <prefix>\n" .
    "\n\tThis requires that kraken2 and kraken2-build be installed\n" .
    "\t\tkraken2 command that would be used: '$kcmd'\n" .
    "\t\tkraken2-build command that would be used: '$kbcmd'\n";
my ($mira_t, $mira_txt, $fasta, $read, $prefix, @rest) = @ARGV;
die "The '-t' flag is missing from the invocation\n$usage" unless ($mira_t && $mira_t eq "-t");
die "The txt from the '-t' flag is missing from the invocation\n$usage" unless ($mira_txt && $mira_txt eq "txt");
die "ERROR: Missing argument: no bait fasta defined\n" . $usage unless $fasta;
die "ERROR: Missing argument: no read file defined\n" . $usage unless $read;
die "ERROR: Missing argument: no prefix defined\n" . $usage unless $prefix;
die "$prefix\.txt already exists in the directory\n$usage" if (-e "$prefix\.txt");
die "There were to many arguments specified at invocation\n$usage" if @rest;

# # Check read file: gziped or not?
# my $read_handle;
# my $gzipped;
# system("gzip -l $read 2>/dev/null >/dev/null");
# if ($? == 0) {
#     $gzipped = "yes";
# }


# Identify reads that match
# Open output file
open(my $out, '>', "$prefix\.txt") || die $!;

# Set up kraken protein DB
my $dir = "prot-kmer-db";

# Create input and directory for the DB
if (-d $dir) {
    `rm -r $dir`;
    print STDERR "Removing old DB directory ($dir)\n";
}
mkdir $dir;
open my $bait, '<', $fasta || die $!;
open my $input, '>', 'input.faa' || die $1;;
while(<$bait>) {
    s/>(\S+)/>$1|kraken:taxid|194397/;
    print {$input} $_;
}

# Set up taxonomy for the DB
mkdir $dir . "/taxonomy";
my $path = $dir . "/taxonomy/";

# Add minimal taxonomy files for Fusarium graminearum dsRNA mycovirus-1
#  The taxonomy information is incorrect, but it is only needed to set up
#  a kraken2 protein DB that will be used for baiting only, not for true classification.
open my $names, '>', $path . 'names.dmp';
print {$names} "",
"1	|	all	|		|	synonym	|\n",
"1	|	root	|		|	scientific name	|\n",
"10239	|	Vira	|		|	synonym	|\n",
"10239	|	Viridae	|		|	synonym	|\n",
"10239	|	viruses	|		|	blast name	|\n",
"10239	|	Viruses	|		|	scientific name	|\n",
"35325	|	dsRNA nonenveloped viruses	|		|	equivalent name	|\n",
"35325	|	dsRNA viruses	|		|	scientific name	|\n",
"39780	|	unclassified dsRNA viruses	|		|	scientific name	|\n",
"194397	|	Fusarium graminearum dsRNA mycovirus-1	|		|	scientific name	|\n",
"2559587	|	Riboviria	|		|	scientific name	|\n",
"2559587	|	RNA viruses and viroids	|		|	common name	|\n",
"2559587	|	RNA viruses	|		|	common name	|\n";
close $names;

open my $nodes, '>', $path . 'nodes.dmp';
print {$nodes} "",
"1	|	1	|	no rank	|		|	8	|	0	|	1	|	0	|	0	|	0	|	0	|	0	|		|\n",
"10239	|	1	|	superkingdom	|		|	9	|	0	|	1	|	0	|	0	|	0	|	0	|	0	|		|\n",
"2559587	|	10239	|	clade	|		|	9	|	1	|	1	|	1	|	0	|	1	|	0	|	0	|		|\n",
"35325	|	2559587	|	clade	|		|	9	|	1	|	1	|	1	|	0	|	1	|	0	|	0	|		|\n",
"39780	|	35325	|	no rank	|		|	9	|	1	|	1	|	1	|	0	|	1	|	1	|	0	|		|\n",
"194397	|	39780	|	species	|	FG	|	9	|	1	|	1	|	1	|	0	|	1	|	1	|	0	|		|\n";

# Commands to build the DB
`$kbcmd --add-to-library input.faa --protein --db $dir`;
`$kbcmd --db $dir --protein --build`;

# Run kraken protein analysis and print output
open(my $run, "$kcmd --db $dir $read |") || die $!;
while(<$run>){
    # Skip reads that did not match the reference/bait
    next if /^U/;
    # collect and print read IDs to output file
    s/\R//g;
    # Only the second column (separator = \t) is needed
    my ($nul, $read_id) = split/\t/;
    print {$out} $read_id . "\n";
}

# Done and exit correctly
exit 0;
