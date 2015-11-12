#!/usr/bin/perl -w
use strict;
use autodie;


# Check which version they want to configure
#  minimal or normal
my $response;
if (@ARGV == 0) {
    print "This is a script to test the installation of prerequisites for the programs in\n" .
	"the GRAbB package, and to configure the scripts before usage.\n";
    print "This script may be run with the following arguments:\n";
    print "\tall\t- configure GRAbB and the helper programs\n";
    print "\tdef\t- configure GRAbB\n";
    print "\tmin\t- configure minmal_GRAbB\n";
    print "\thelp\t- configure the helper programs\n\n";
    print "Which option would you like to run? ";
    $response = <>;
    push @ARGV, $response;
}

my ($normal, $min, $help);
for (@ARGV) {
    if (/all/) {
	$normal++;
	$help++;
    } elsif (/help/) {
	$help++;
    } elsif (/min/) {
	$min++;
    } elsif (/def/) {
	$normal++;
    }
}

print "Configuring GRAbB package\n";

# Variables that are true if given element is functional
my ($bait, $collect, $assembler, $exonerate, $prinseq);

#my $prog; # store the name of the current program that is being checked
# Check baiting program
my $call ='-t txt for_testing/assembly.fas for_testing/read1.fastq temp 2>err >out';
if ($normal) {
    # Check mirabait
    my $prog = 'mirabait';
    print "\tTesting $prog\n";
    
    my $cmd;
    # Test if the program is available through the path
    if (`which $prog`) {
	print "\t\t$prog is found in the path\n";
	$cmd = $prog;

	# Check if it is working correctly
	system("$cmd $call");
	if ($? == 0 && -s 'temp.txt') {
	    print "\t\t$prog (found in the path) is working correctly\n";
	    # Save the command that works
	    $bait = $cmd;
	} else {
	    print "\t\t$prog (found in the path) is not working correctly\n";
	}
	# Remove possible output files
	for ('out', 'err', 'temp.txt', 'hashstat.bin') {
	    unlink $_ if -e $_;
	}
    } else {
	print "\t\t$prog is not found in the path\n";
    }
}
if ($normal && not $bait) {
    # Check the mirabait included in the package
    my $prog = 'mirabait';
    if (-x '3rd_party_programs/mira_4.0rc4_linux-gnu_x86_64_static/bin/mirabait') {
	print "\t\t$prog executable found inside the package\n";
	# Use the absolute path of the executable
	my $cmd = $ENV{"PWD"} . "/" . '3rd_party_programs/mira_4.0rc4_linux-gnu_x86_64_static/bin/mirabait';

	# Check if it is working correctly
	system("$cmd $call");
	if ($? == 0 && -s 'temp.txt') {
	    print "\t\t$prog (found in the package) is working correctly\n";
	    # Save the command that works
	    $bait = $cmd;
	} else {
	    print "\t\t$prog (found in the package) is not working correctly\n";
	}
	# Remove possible output files
	for ('out', 'err', 'temp.txt', 'hashstat.bin') {
	    unlink $_ if -e $_;
	}
    } else {
	print "Could not locate executable for $prog inside the package\n";
    }
}
if (not $bait) {
    print "\tConfiguration script failed to locate a working version of mirabait\n" if $normal;
    my $prog = 'kmer_bait.pl';

    print "\tTesting $prog\n";
    
    my $cmd;
    # Test if the program is available through the path
    if (`which $prog`) {
	print "\t\t$prog is found in the path\n";
	$cmd = $prog;

	# Check if it is working correctly
	system("$cmd $call");
	if ($? == 0 && -s 'temp.txt') {
	    print "\t\t$prog (found in the path) is working correctly\n";
	    # Save the command that works
	    $bait = $cmd;
	} else {
	    print "\t\t$prog (found in the path) is not working correctly\n";
	}
	# Remove possible output files
	for ('out', 'err', 'temp.txt', 'hashstat.bin') {
	    unlink $_ if -e $_;
	}
    } else {
	print "\t\t$prog is not found in the path\n";
    }
    # Try the file in the package
    if (not $bait) {
	if (-x 'perl_programs/kmer_bait.pl') {
	    print "\t\t$prog executable found inside the package\n";
	    # Use the absolute path of the executable
	    my $cmd = $ENV{"PWD"} . "/" . 'perl_programs/kmer_bait.pl';
	    
	    # Check if it is working correctly
	    system("$cmd $call");
	    if ($? == 0) {
		print "\t\t$prog (found in the package) is working correctly\n";
		# Save the command that works
		$bait = $cmd;
	    } else {
		print "\t\t$prog (found in the package) is not working correctly\n";
	    }
	    # Remove possible output files
	    for ('out', 'err', 'temp.txt') {
		unlink $_ if -e $_;
	    }
	} else {
	    print "\t\tCould not locate executable for $prog inside the package\n";
	    if (-e 'perl_programs/kmer_bait.pl') {
		print "\t\t$prog source found inside the package\n";
		# Use the absolute path of the executable
		my $cmd = "perl " . $ENV{"PWD"} . "/" . 'perl_programs/kmer_bait.pl';

		# Check if it is working correctly
		system("$cmd $call");
		if ($? == 0 && -s 'temp.txt') {
		    print "\t\t$prog (found in the package) is working correctly\n";
		    # Save the command that works
		    $bait = $cmd;
		} else {
		    print "\t\t$prog (found in the package) is not working correctly\n";
		}
		# Remove possible output files
		for ('out', 'err', 'temp.txt') {
		    unlink $_ if -e $_;
		}
	    } else {
		print "\t\t$prog source could not be located\n";
	    }
	}
    }
}

# Check collecting program
$call = 'for_testing/read1.fastq for_testing/found.txt 2>err >out';
if ($normal) {
    # Check seqtk
    my $prog = 'seqtk';
    print "\tTesting $prog\n";
    
    my $cmd;
    # Test if the program is available through the path
    if (`which $prog`) {
	print "\t\t$prog is found in the path\n";
	$cmd = $prog;

	# Check if it is working correctly
	system("$cmd subseq $call");
	if ($? == 0 && -s 'out') {
	    print "\t\t$prog (found in the path) is working correctly\n";
	    # Save the command that works
	    $collect = $cmd;
	} else {
	    print "\t\t$prog (found in the path) is not working correctly\n";
	}
	# Remove possible output files
	for ('out', 'err') {
	    unlink $_ if -e $_;
	}
    } else {
	print "\t\t$prog is not found in the path\n";
    }
}
if (not $collect) {
    print "\tConfiguration script failed to locate a working version of seqtk\n" if $normal;
    my $prog = 'create_readpool.pl';

    print "\tTesting $prog\n";
    
    my $cmd;
    # Test if the program is available through the path
    if (`which $prog`) {
	print "\t\t$prog is found in the path\n";
	$cmd = $prog;

	# Check if it is working correctly
	system("$cmd $call");
	if ($? == 0 && -s 'out') {
	    print "\t\t$prog (found in the path) is working correctly\n";
	    # Save the command that works
	    $collect = $cmd;
	} else {
	    print "\t\t$prog (found in the path) is not working correctly\n";
	}
	# Remove possible output files
	for ('out', 'err') {
	    unlink $_ if -e $_;
	}
    } else {
	print "\t\t$prog is not found in the path\n";
    }
    # Try the file in the package
    if (not $collect) {
	if (-x "perl_programs/$prog") {
	    print "\t\t$prog executable found inside the package\n";
	    # Use the absolute path of the executable
	    my $cmd = $ENV{"PWD"} . "/" . "perl_programs/$prog";
	    
	    # Check if it is working correctly
	    system("$cmd $call");
	    if ($? == 0) {
		print "\t\t$prog (found in the package) is working correctly\n";
		# Save the command that works
		$collect = $cmd;
	    } else {
		print "\t\t$prog (found in the package) is not working correctly\n";
	    }
	    # Remove possible output files
	    for ('out', 'err') {
		unlink $_ if -e $_;
	    }
	} else {
	    print "\t\tCould not locate executable for $prog inside the package\n";
	    if (-e "perl_programs/$prog") {
		print "\t\t$prog source found inside the package\n";
		# Use the absolute path of the executable
		my $cmd = "perl " . $ENV{"PWD"} . "/" . "perl_programs/$prog";

		# Check if it is working correctly
		system("$cmd $call");
		if ($? == 0 && -s 'out') {
		    print "\t\t$prog (found in the package) is working correctly\n";
		    # Save the command that works
		    $collect = $cmd;
		} else {
		    print "\t\t$prog (found in the package) is not working correctly\n";
		}
		# Remove possible output files
		for ('out', 'err', 'temp.txt', 'hashstat.bin') {
		    unlink $_ if -e $_;
		}
	    } else {
		print "\t\t$prog source could not be located\n";
	    }
	}
    }
}

# Check assembler program
# Check Edena
# Check Velvet
# Check Exonerate






# For PRINSEQ-lite the following modules are needed
if ($help) {
    # Test each module that is needed for PRINSEQ-lite
    # Collect the missing modules at the end
    my @missing;
    if (eval {require Getopt::Long; 1;} ne 1) {
	push @missing, "Getopt::Long";
    }
    if (eval {require Pod::Usage; 1;} ne 1) {
	push @missing, "Pod::Usage";
    }
    if (eval {require File::Temp; 1;} ne 1) {
	push @missing, "File::Temp";
    }
    if (eval {require Fcntl; 1;} ne 1) {
	push @missing, "Fcntleval";
    }
    if (eval {require Digest::MD5; 1;} ne 1) {
	push @missing, "Digest::MD5";
    }
    if (eval {require Cwd; 1;} ne 1) {
	push @missing, "Cwd";
    }
    if (eval {require List::Util; 1;} ne 1) {
	push @missing, "List::Util";
    }
    
    if (@missing) {
	print "The following Perl modules are not installed, but are required for PRINSEQ-lite:\n\t\t";
	print join("\n\t\t", @missing);
	print "\n\tUse cpan or cpanm to install these modules (http://www.cpan.org/modules/INSTALL.html)\n";
    } else {
	# Configure helper scripts
	print "All the Perl modules are installed that are needed to run PRINSEQ-lite.\n";
    }
}
