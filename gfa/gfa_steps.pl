#!/usr/bin/perl -w
use strict;


my $ref = $ARGV[0];
my $gfa = $ARGV[1];
my $result = $ARGV[2];
my $fasta;
$fasta = $ARGV[3] if scalar @ARGV > 3;

# Create fasta file
my $target = "gfa.fas";
`gfa2fasta.pl $gfa >$target`;

my $ryostring = '--ryo "ryo\t%ti\t%tl\t%et\t%ei\n"';
my @result = `exonerate $ref $target $ryostring`;

@result = map{s/^ryo\t//; $_} grep {/^ryo\t/} @result;

my %hit;
for (@result) {
    my ($id, $len, $align, $ident) = split/\t/;
    my $sim = $ident / $align;
    if ($ident > 31) {
	$hit{$id}++;
#	print
    }
}

my $get = join("\t", sort keys %hit);

print $get . "\n";

`gfa_extract.pl $gfa $get >$result`;
`gfa2fasta.pl $result >$fasta` if $fasta;
