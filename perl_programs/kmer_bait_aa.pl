#!/usr/bin/perl -w
use strict;
use autodie;

# This is a perl alternative for using mirabait to bait reads using exact 31-mer matching
# Invocation by GRAbB.pl:
#    $bait_cmd -t txt $$bait_ref $read $mira_temp\_$i 2>&1 >>mirabait.log

# Input:
#        "-t"
#        "txt"
#        Reference (Fasta)
#        Read file (Fasta or Fastq)
#        prefix (prefix for the output file: prefix.txt)

    my %codons = (
	'TCA'=>'S', #Serine
	'TCC'=>'S', #Serine
	'TCG'=>'S',  #Serine
	'TCT'=>'S', #Serine
	'TTC'=>'F', #Phenylalanine
	'TTT'=>'F', #Phenylalanine
	'TTA'=>'L', #Leucine
	'TTG'=>'L', #Leucine
	'TAC'=>'Y', #Tyrosine
	'TAT'=>'Y', #Tyrosine
	'TAA'=>'_', #Stop
	'TAG'=>'_', #Stop
	'TGC'=>'C', #Cysteine
	'TGT'=>'C', #Cysteine
	'TGA'=>'_', #Stop
	'TGG'=>'W', #Tryptophan
	'CTA'=>'L', #Leucine
	'CTC'=>'L', #Leucine
	'CTG'=>'L', #Leucine
	'CTT'=>'L', #Leucine
	'CCA'=>'P', #Proline
	'CAT'=>'H', #Histidine
	'CAA'=>'Q', #Glutamine
	'CAG'=>'Q', #Glutamine
	'CGA'=>'R', #Arginine
	'CGC'=>'R', #Arginine
	'CGG'=>'R', #Arginine
	'CGT'=>'R', #Arginine
	'ATA'=>'I', #Isoleucine
	'ATC'=>'I', #Isoleucine
	'ATT'=>'I', #Isoleucine
	'ATG'=>'M', #Methionine
	'ACA'=>'T', #Threonine
	'ACC'=>'T', #Threonine
	'ACG'=>'T', #Threonine
	'ACT'=>'T', #Threonine
	'AAC'=>'N', #Asparagine
	'AAT'=>'N', #Asparagine
	'AAA'=>'K', #Lysine
	'AAG'=>'K', #Lysine
	'AGC'=>'S', #Serine
	'AGT'=>'S', #Serine
	'AGA'=>'R', #Arginine
	'AGG'=>'R', #Arginine
	'CCC'=>'P', #Proline
	'CCG'=>'P', #Proline
	'CCT'=>'P', #Proline
	'CAC'=>'H', #Histidine
	'GTA'=>'V', #Valine
	'GTC'=>'V', #Valine
	'GTG'=>'V', #Valine
	'GTT'=>'V', #Valine
	'GCA'=>'A', #Alanine
	'GCC'=>'A', #Alanine
	'GCG'=>'A', #Alanine
	'GCT'=>'A', #Alanine
	'GAC'=>'D', #Aspartic Acid
	'GAT'=>'D', #Aspartic Acid
	'GAA'=>'E', #Glutamic Acid
	'GAG'=>'E', #Glutamic Acid
	'GGA'=>'G', #Glycine
	'GGC'=>'G', #Glycine
	'GGG'=>'G', #Glycine
	'GGT'=>'G', #Glycine
    );


# Test correct invocation
my $usage = "Usage\n\t$0 -t txt <reference.fas> <read.fastq|read.fasta> <prefix>\n";
my ($mira_t, $mira_txt, $fasta, $read, $prefix, @rest) = @ARGV;
die "The '-t' flag is missing from the invocation\n$usage" unless $mira_t eq "-t";
die "The txt from the '-t' flag is missing from the invocation\n$usage" unless $mira_txt eq "txt";
die "$prefix\.txt already exists in the directory\n$usage" if (-e "$prefix\.txt");
die $usage unless $fasta;
die $usage unless $read;
die "There were to many arguments specified at invocation\n$usage" if @rest;

# read fasta
my %fas_data;
my @ids;
&read_fasta(\%fas_data, \@ids, $fasta);

# Generate all the k-mers for all the sequences in the fasta file
my $k_length = 10;
my %kmer;
for (@ids) {
    &gen_kmer($fas_data{$_}, \%kmer, $k_length);
}

# Identify reads that match
# Open output file
open(my $out, '>', "$prefix\.txt") || die $!;

# open read file (gziped or not)
my $read_handle;
#system("zcat $read | head -1 >/dev/null");
system("gzip -l $read 2>/dev/null >/dev/null");
if ($? == 0) {
    open($read_handle, sprintf("zcat %s |", $read)) || die $!;
} else {
    open($read_handle, '<', $read) || die $!;
}

my $read_id;
my $read_seq;
# is it fasta or fastq
# Get the first line
my $first = <$read_handle>;
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
	    /\@(\S+)/;
	    $read_id = $1;
	    $read_seq = "";
	} elsif ($n%4 == 2) {
	    # get read sequence
	    s/\R//g;
	    $read_seq = $_;
	    # check for kmer match
	    print {$out} "$read_id\n" if &aa_kmer_match($read_seq, \%kmer, $k_length);
	}
    }
} else {
    for (<$read_handle>) {
	# save read id
	if (/^>(\S+)/) {
	    my $new_id = $1;
	    print {$out} "$read_id\n" if &aa_kmer_match($read_seq, \%kmer, $k_length);
	    $read_id = $new_id;
	    $read_seq = "";
	} else {
	    # get read sequence
	    s/\R//g;
	    $read_seq .= $_;
	}
    }
    print {$out} "$read_id\n" if &aa_kmer_match($read_seq, \%kmer, $k_length);
}


exit 0;

#=======Subroutines====================================

sub read_fasta {
    # This loads the sequences into a hash and an array
    my ($hash, $list, $file) = @_;
    my $in;
    if ($file) {
	open $in, '<', $file || die $!;
    } else {
	$in = *STDIN;
    }
    my $seq_id;
    for (<$in>) {
	# Skip empty lines
	next if /^\s*$/;
	# Check wheter it is an id line
	if (/>(.*)/) {
	    # Save the id and the definition and store it in the array
	    # The sequence identifier is not needed
	    $seq_id++;
	    push @$list, $seq_id;
	} else {
	    # If there was no id lines before this then throw an error
	    unless (defined $seq_id) {
		print "Format error! Check the file!\n";
		last;
	    }
	    # Remove lineendings and white space
	    s/\R//g;
	    s/\s+//g;
	    # Add to the sequence
	    $hash->{$seq_id} .= $_;
	}
    }
    close $in;
}


sub gen_kmer {
# generates kmers from a sequence
    my $seq = shift;
    my $kmer_ref = shift;
    my $k = shift;
    my $len = length $seq;
    for (0..($len-$k)) {
	my $kmer = substr $seq, $_, $k;
	$$kmer_ref{$kmer}++;
#	$$kmer_ref{&reverse_seq($kmer)}++; # Not needed for proteins, no complement
    }
}

sub aa_kmer_match {
#  &kmer_match($read_seq, \%kmer, $k_length)
   my ($seq, $kmer_ref, $k) = @_;
   my $aa = &translate($seq);
   return 1 if &kmer_match($aa, $kmer_ref, $k);
   $aa = &translate(substr $seq, 1);
   return 1 if &kmer_match($aa, $kmer_ref, $k);
   $aa = &translate(substr $seq, 2);
   return 1 if &kmer_match($aa, $kmer_ref, $k);
   $seq = &reverse_seq($seq);
   $aa = &translate($seq);
   return 1 if &kmer_match($aa, $kmer_ref, $k);
   $aa = &translate(substr $seq, 1);
   return 1 if &kmer_match($aa, $kmer_ref, $k);
   $aa = &translate(substr $seq, 2);
   return 1 if &kmer_match($aa, $kmer_ref, $k);
   return 0;
}

sub kmer_match {
#  &kmer_match($read_seq, \%kmer, $k_length)
   my ($seq, $kmer_ref, $k) = @_;
   my $len = length $seq;
   for (0..($len-$k)) {
       my $kmer = substr $seq, $_, $k;
       return 1 if $$kmer_ref{$kmer};
   }
   # Return false since no k-mers were matched
   return 0;
}

sub reverse_seq {
    # Reverse complements the sequences
    my ($seq) = @_;
    # Reverse the sequnce
    my $complement = reverse $seq;
    # Complement the sequence
    $complement =~ tr/ACGTacgtWwMmRrSsKkYyBbVvDdHh/TGCAtgcaWwKkYySsMmRrVvBbHhDd/;
    return $complement;
}


sub translate{
    # Takes string as input, which is the dna sequence
    my ($seq) = @_;
    # Split the sequence to triplets
    my @triplet = $seq =~ /.../g;
    my $prot = "";
    for (@triplet) {
	# Add X if IUPAC base is found in the codon
        $prot .= $codons{$_} ? $codons{$_} : "X";
    }
    # Print the sequences to STDOUT
    return $prot;
}
