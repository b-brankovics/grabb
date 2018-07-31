#!/usr/bin/perl -w

use strict;
#use Data::Dumper;


# Store info on 
my %map;
my %gfa;


# Turn on for debugging
my $debug = 0;


my $min = 0;
# Read the input
if (@ARGV) {
    for (@ARGV) {
	if (/^-?-m(?:in)?=(\d+)$/) {
	    $min = $1;
	    &read_gfa(\%map, \%gfa) if scalar(@ARGV) == 1;
	} else {
	    &read_gfa(\%map, \%gfa, $_);
	}
    }
} else {
    	&read_gfa(\%map, \%gfa);
}


# Create fasta data for segments with sufficient lengths
my @ids;
my %fas_data;    
&gfa2fasta(\%gfa, \@ids, \%fas_data, $min);

# Print the FASTA sequences to STDOUT
&print_fasta(\%fas_data, \@ids);




# Subrutines

sub read_gfa {
    my ($map, $gfa, $file) = @_;
    # Use STDIN if file is '-'
    $file = undef if $file && $file eq '-';
    my $in;
    if ($file && -e $file) {
	open $in, '<', $file || die $!;
    } else {
	$in = *STDIN;
    }

    for (<$in>) {
	# Remove line endings
	s/\R//g;
	
	# split the line
	my @tab = split/\t/;
	my $type = $tab[0];
	
	if ($type eq "S") {
	    # Segement entry
	    print "" . join("\t", $tab[0], $tab[1]) . "\n" if ($debug);
	    my $id = $tab[1];
	    # If already has a mapping
	    #  - then append to appropriate %gfa entry
	    #  ELSE create new %gfa entry with current line + create mapping
	    if ($map->{$id}) {
		print "append to gfa (" . $map->{$id} . ")\n" if $debug;
		push @{ $gfa->{ $map->{$id} } }, $_;
	    } else {
		print "create new gfa ($id)\n" if $debug;
		$map->{$id} = $id;
		$gfa->{$id} = [$_];
	    }
	    
	} elsif ($type eq "L") {
	    # Link entry
	    print "" . join("\t", @tab[0..4]) . "\n" if ($debug);
	    my $a = $tab[1];
	    my $b = $tab[3];
	    my $id;
	    if ($map->{$a}) {
		$id = $map->{$a};
		if ($map->{$b}) {
		    # Check if they are already the same
		    unless ($id eq $map->{$b}) {
			print "$a => " . $map->{$a} . "\n$b => " . $map->{$b} ."\n" if $debug;
			# merge the two arrays, delete $map{$b} in %gfa, update map
			push @{ $gfa->{$id} }, @{ $gfa->{ $map->{$b} } };
			delete $gfa->{ $map->{$b} };
			&update_map($map, $map->{$b}, $id);
			$map->{$b} = $id;
		    }
		} else {
		    # Create a "link"
		    $map->{$b} = $id;
		}
	    } else {
		if ($map->{$b}) {
		    $id = $map->{$b};
		    # Create a "link"
		    $map->{$a} = $id;
		} else {
		    # Create new entries
		    $id = $a;
		    $map->{$a} = $id;
		    $map->{$b} = $id;
		    $gfa->{$id} = [];
		}
	    }
	    # Append line
	    push @{ $gfa->{$id} }, $_;
	}
	
    }
}

sub update_map {
    # subrutine for updating %map
    # When removing an element from %gfa, we need to update the old pointers from $old to $new
    my ($hash, $old, $new) = @_;
    for (keys %$hash) {
	$hash->{$_} = $new if $hash->{$_} eq $old;
    }    
}

sub gfa2fasta {
    # Create fasta info from gfa (segments)
    my ($gfa, $ids, $fas_data, $min) = @_;
    $min = 0 unless $min;
    for my $graphid (keys %$gfa) {
	# Go through the graph lines
	for (@{ $gfa->{$graphid} }) {
	    my ($type, $id, $seq) = split/\t/;
	    # Skip if link or too short
	    if ($type eq "S" && length($seq) >= $min) {
		push @$ids, $id;
		$fas_data->{$id} = $seq;
	    }
	}
    }
}

sub to_fasta {
    # Return a fasta formated string
    my ($seq_name, $seq, $len) = @_;
    # default to 60 characters of sequence per line
    $len = 60 unless $len;
    # Print ID line
    my $formatted_seq = ">$seq_name\n";
    # Add sequence lines with $len length
    while (my $chunk = substr($seq, 0, $len, "")) {
	$formatted_seq .= "$chunk\n";
    }
    return $formatted_seq;
}

sub print_fasta {
    # Print all the sequences to STDOUT in FASTA format
    my ($hash, $list) = @_;
    for (@$list) {
	print &to_fasta($_, $hash->{$_});
    }
}
