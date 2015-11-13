#!/usr/bin/perl -w
use strict;
use autodie;

# This is a perl alternative for using seqtk to collect reads that were matched by kmer matching
# Invocation by GRAbB.pl:
#    $collect_cmd $$read_ref $$list_ref >$$readpool_ref

# Input:
#        Read file (Fasta or Fastq)
#        List file (txt)

# Test correct invocation
my $usage = "Usage\n\t$0 <read.fastq|read.fasta> <list>\n";
my ($readf, $list, @rest) = @ARGV;
die $usage unless $readf;
die $usage unless $list;
die "There were to many arguments specified at invocation\n$usage" if @rest;

# Get list of identifiers
my %ids;
open(my $txt, '<', $list) || die $!;
for (<$txt>) {
    s/\R//g;
    $ids{$_}++;
}

# Print reads that have matched

# open read file (gziped or not)
my $read_handle;
#system("zcat $read | head -1 >/dev/null");
system("gzip -l $readf 2>/dev/null >/dev/null");
if ($? == 0) {
    open($read_handle, sprintf("zcat %s |", $readf)) || die $!;
} else {
    open($read_handle, '<', $readf) || die $!;
}

my $read_id;
my $read;
# is it fasta or fastq
# Get the first line
my $first = <$read_handle>;
$read = $first;
my $fastq;
if ($first =~ /^>/) {
    # It is fasta
    $fastq = undef;
} elsif ($first =~ /^\@/) {
    # It is fastq
    $fastq++;
} else {
    die "First line of the read file is neither FASTA nor FASTQ format\n";
}

# Get the first read id
$first =~ /^[>\@](\S+)/;
$read_id = $1;
if ($fastq) {
    # Keep track of line numbers
    my $n = 1;
    for (<$read_handle>) {
	$n++;
	# save read id
	if ($n%4 == 1) {
	    print $read if $ids{$read_id};
	    /\@(\S+)/;
	    $read_id = $1;
	    $read = $_;
	} else {
	    # get read
	    $read .= $_;
	}
    }
    print $read if $ids{$read_id};
} else {
    for (<$read_handle>) {
	# save read id
	if (/^>(\S+)/) {
	    my $new_id = $1;
	    print $read if $ids{$read_id};
	    $read_id = $new_id;
	    $read = $_;
	} else {
	    # get read
	    $read .= $_;
	}
    }
    print $read if $ids{$read_id};
}


exit 0;
