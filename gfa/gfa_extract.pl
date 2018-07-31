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

my @requests; # Collect the regular expresions for extraction

if (@ARGV) {
    for (@ARGV) {
	if (-f $_ || $_ eq "-") {
	    &read_gfa(\%map, \%gfa, $_);
	} else {
	    push @requests, $_;
	}
    }
}
unless ( scalar(keys %gfa) ) {
    &read_gfa(\%map, \%gfa);
}


my %need;
for (@requests) {
    $need{ $map{$_} }++;
}

for (keys %gfa) {
    if ($need{$_}) {
	my $ref = $gfa{$_};
	for (@$ref) {
	    print "$_\n";
	}
    }
}



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
