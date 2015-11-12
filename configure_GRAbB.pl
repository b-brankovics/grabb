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
my ($bait, $collect, $edena, $velveth, $velvetg, $exonerate, $prinseq);

{# Check baiting program
    my @input = ('for_testing/assembly.fas', 'for_testing/read1.fastq');
    my $call ="-t txt @input temp 2>err >out";
    my $prog = 'mirabait';
    my $out = 'temp.txt';
    my @trash = ('out', 'err', 'temp.txt', 'hashstat.bin');
    my $path = '3rd_party_programs/mira_4.0rc4_linux-gnu_x86_64_static/bin/mirabait';

    my $working;

    my $good_to_go = 0;
    for (@input) {
	$good_to_go++ if -s $_;
    }

    unless ($good_to_go == scalar @input) {
	print "\tCould not locate test files for $prog\n\t\tMake sure that configure script is run inside the package folder\n";
	last;
    }

    $working = &test_path($prog, $call, $out, @trash) if $normal;
    $working = &test_exe($prog, $path, $call, $out, @trash) if ($normal && not $working);

    print "\tConfiguration script failed to locate a working version of $prog\n" if ($normal && not $working);
    $prog = 'kmer_bait.pl' unless $working;
    $working = &test_path($prog, $call, $out, @trash) unless $working;

    # Try the file in the package
    $path = "perl_programs/$prog" unless $working;
    $working =  &test_exe($prog, $path, $call, $out, @trash) unless $working;
    $working = &test_perl($prog, $path, $call, $out, @trash) unless $working;

    $bait = $working;
}

{# Check collecting program
    my @input = ('for_testing/read1.fastq', 'for_testing/found.txt');
    my $call ="@input 2>err >out";
    my $prog = 'seqtk subseq';
    my $out = 'out';
    my @trash = ('out', 'err');
    my $path = '';
    my $working;

    my $good_to_go = 0;
    for (@input) {
	$good_to_go++ if -s $_;
    }

    unless ($good_to_go == scalar @input) {
	print "\tCould not locate test files for $prog\n\t\tMake sure that configure script is run inside the package folder\n";
	last;
    }

    $working = &test_path($prog, $call, $out, @trash) if $normal;

    print "\tConfiguration script failed to locate a working version of $prog\n" if ($normal && not $working);
    $prog = 'create_readpool.pl' unless $working;
    $working = &test_path($prog, $call, $out, @trash) unless $working;

    # Try the file in the package
    $path = "perl_programs/$prog" unless $working;
    $working =  &test_exe($prog, $path, $call, $out, @trash) unless $working;
    $working = &test_perl($prog, $path, $call, $out, @trash) unless $working;

    $collect = $working;
}

{# Check Edena
    my @input = ('for_testing/read1.fastq', 'for_testing/read2.fastq');
    my $call1 ="-paired @input 2>err >out";
    my $call2 ="-e out.ovl 2>err >out";
    my $prog = 'edena';
    my $out = 'out.ovl';
    my $out2 = 'out_contigs.fasta';
    my @trash = ('out', 'err', 'out_overlapping.log','out.ovl');
    my $path = '3rd_party_programs/EdenaV3.131028/bin/edena';

    my $working;

    my $good_to_go = 0;
    for (@input) {
	$good_to_go++ if -s $_;
    }

    unless ($good_to_go == scalar @input) {
	print "\tCould not locate test files for $prog\n\t\tMake sure that configure script is run inside the package folder\n";
	last;
    }

    print "\tTesting $prog\n";
    # Test if the program is available through the path
    my $print = "\t\t$prog (found in the path)";
    if (`which $prog`) {
	print "\t\t$prog is found in the path\n";
	my $cmd = $prog;

	# Check if it is working correctly
	# Command, arguments, output, text, to be removed
	system("$cmd $call1");
	if ($? == 0 && -s $out) {
	    print "$print step 1 is working correctly\n";
	    # Save the command that works
	    system("$cmd $call2");
	    if ($? == 0 && -s $out2) {
		print "$print step 2 is working correctly\n";
		# Save the command that works
		$working = $cmd;
	    } else {
		print "$print step 2 is not working correctly\n";
	    }
	    for (glob("out*")) {
		unlink $_ if -e $_;
	    }
	} else {
	    print "$print step 1 is not working correctly\n";
	}
	# Remove possible output files
	for (@trash) {
	    unlink $_ if -e $_;
	}
    } else {
	print "\t\t$prog is not found in the path\n";
    }

    unless ($working) {    
	
	$print = "\t\t$prog (found in the package)";
	if (-x $path) {
	    print "\t\t$prog is found in the package\n";
	    my $cmd = $ENV{"PWD"} . '/' . $path;

	    # Check if it is working correctly
	    system("$cmd $call1");
	    if ($? == 0 && -s $out) {
		print "$print step 1 is working correctly\n";

		system("$cmd $call2");
		if ($? == 0 && -s $out2) {
		    print "$print step 2 is working correctly\n";
		    # Save the command that works
		    $working = $cmd;
		} else {
		    print "$print step 2 is not working correctly\n";
		}
		for (glob("out*")) {
		    unlink $_ if -e $_;
		}
	    } else {
		print "$print step 1 is not working correctly\n";
	    }
	    # Remove possible output files
	    for (@trash) {
		unlink $_ if -e $_;
	    }
	} else {
	    print "\t\t$prog is not found inside the package\n";
	}
    }

    print "\tConfiguration script failed to locate a working version of $prog\n" unless $working;

    $edena = $working;
}

{# Check Velveth
    my @input = ('for_testing/read1.fastq', 'for_testing/read2.fastq');
    my $call ="velvet 31 -shortPaired -fastq -separate @input 2>err >out";
    my $prog = 'velveth';
    my $out = 'velvet/Roadmaps';
    my $path = '3rd_party_programs/velvet_1.2.10/velveth';

    my $working;

    my $good_to_go = 0;
    for (@input) {
	$good_to_go++ if -s $_;
    }

    unless ($good_to_go == scalar @input) {
	print "\tCould not locate test files for $prog\n\t\tMake sure that configure script is run inside the package folder\n";
	last;
    }

    print "\tTesting $prog\n";
    # Test if the program is available through the path
    my $print = "\t\t$prog (found in the path)";
    if (`which $prog`) {
	print "\t\t$prog is found in the path\n";
	my $cmd = $prog;

	# Check if it is working correctly
	system("$cmd $call");
	if ($? == 0 && -s $out) {
	    print "$print is working correctly\n";
	    # Save the command that works
	    $working = $cmd;
	} else {
	    print "$print is not working correctly\n";
	}
    } else {
	print "\t\t$prog is not found in the path\n";
    }

    unless ($working) {    
	
	$print = "\t\t$prog (found in the package)";
	if (-x $path) {
	    print "\t\t$prog is found in the package\n";
	    my $cmd = $ENV{"PWD"} . '/' . $path;

	    # Check if it is working correctly
	    system("$cmd $call");
	    if ($? == 0 && -s $out) {
		print "$print is working correctly\n";
		# Save the command that works
		$working = $cmd;
	    } else {
		print "$print step 1 is not working correctly\n";
	    }
	} else {
	    print "\t\t$prog is not found inside the package\n";
	}
    }

    print "\tConfiguration script failed to locate a working version of $prog\n" unless $working;

    $velveth = $working;
}

{# Check Velvetg
    my $prog = 'velvetg';
    my $out = 'velvet/contigs.fa';
    my $path = '3rd_party_programs/velvet_1.2.10/velvetg';
    my $call ="velvet >out 2>err";

    my $working;

    unless (-d "velvet") {
	print "\tCould not locate directory produced by velveth\n";
	last;
    }

    print "\tTesting $prog\n";
    # Test if the program is available through the path
    my $print = "\t\t$prog (found in the path)";
    if (`which $prog`) {
	print "\t\t$prog is found in the path\n";
	my $cmd = $prog;

	# Check if it is working correctly
	system("$cmd $call");
	if ($? == 0 && -s $out) {
	    print "$print is working correctly\n";
	    # Save the command that works
	    $working = $cmd;
	} else {
	    print "$print is not working correctly\n";
	}
    } else {
	print "\t\t$prog is not found in the path\n";
    }

    unless ($working) {    
	
	$print = "\t\t$prog (found in the package)";
	if (-x $path) {
	    print "\t\t$prog is found in the package\n";
	    my $cmd = $ENV{"PWD"} . '/' . $path;

	    # Check if it is working correctly
	    system("$cmd $call");
	    if ($? == 0 && -s $out) {
		print "$print is working correctly\n";
		# Save the command that works
		$working = $cmd;
	    } else {
		print "$print step 1 is not working correctly\n";
	    }
	} else {
	    print "\t\t$prog is not found inside the package\n";
	}
    }


    # Remove possible output files
    unlink(glob("velvet/*"));
    rmdir("velvet");

    print "\tConfiguration script failed to locate a working version of $prog\n" unless $working;

    $velvetg = $working;
}

{# Check exonerate
    my @input = ('for_testing/query.fas', 'for_testing/assembly.fas');
    my $call = "@input --percent 5  --model affine:overlap -E --ryo ";
    $call .= '"(%qi %ql %qab %qae %qal %qS) (%ti %tl %tab %tae %tal %tS) (%V) (%et %es %em %pi)\n" ';
    $call .= "--showalignment no --showvulgar no 2>err >out";
    my $prog = 'exonerate';
    my $out = 'out';
    my @trash = ('out', 'err');

    my $working;

    my $good_to_go = 0;
    for (@input) {
	$good_to_go++ if -s $_;
    }

    unless ($good_to_go == scalar @input) {
	print "\tCould not locate test files for $prog\n\t\tMake sure that configure script is run inside the package folder\n";
	last;
    }

    $working = &test_path($prog, $call, $out, @trash) if $normal;

    print "\tConfiguration script failed to locate a working version of $prog\n" if ($normal && not $working);

    $exonerate = $working;
}

# For PRINSEQ-lite the following modules are needed
if ($help) {
    print "Configuring helper_programs\n";
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
	push @missing, "Fcntl";
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
	print "\tThe following Perl modules are not installed, but are required for PRINSEQ-lite:\n\t\t";
	print join("\n\t\t", @missing);
	print "\n\tUse cpan or cpanm to install these modules (http://www.cpan.org/modules/INSTALL.html)\n";
    } else {
	# Configure helper scripts
	print "\tAll the Perl modules are installed that are needed to run PRINSEQ-lite.\n";
	my $working;
	my $cmd = "prinseq-lite.pl";

	if (`which $cmd`) {
	    $working = $cmd;
	    print "\t\tprinseq-lite found in the path\n";
	} else {
	    print "$cmd is not in the path\n";
	    $cmd = $ENV{"PWD"} . '/' . '3rd_party_programs/prinseq-lite-0.20.4/prinseq-lite.pl';
	    if (-x $cmd) {
		print "\t\tprinseq-lite executable is found in the package\n";
		$working = $cmd;
	    } else {
		print "\t\tprinseq-lite executable is not found in the package\n";
		$working = "perl $cmd" if -s $cmd;
		print "\t\tprinseq-lite source is found in the package\n" if $working;
	    }
	}
	$prinseq = $working;
	print "\t\tprinseq-lite source is not found in the package\n" unless $working;
    }
}

# Add all the collected commands to the corresponding scripts

sub test_path {
    my ($prog, $args, $out, @trash) = @_;
    my $final;
    print "\tTesting $prog\n";
    # Test if the program is available through the path
    if (`which $prog`) {
	print "\t\t$prog is found in the path\n";
	my $cmd = $prog;

	# Check if it is working correctly
	$final = &test_cmd($cmd, $args, $out, "\t\t$prog (found in the path)", @trash);
    } else {
	print "\t\t$prog is not found in the path\n";
    }
    return $final;
}

sub test_exe {
    my ($prog, $path, $args, $out, @trash) = @_;
    my $final;
    if (-x $path) {
	print "\t\t$prog executable found inside the package\n";
	# Use the absolute path of the executable
	my $cmd = $ENV{"PWD"} . "/" . $path;

	# Check if it is working correctly
	$final = &test_cmd($cmd, $args, $out, "\t\t$prog (found in the package)", @trash);
    } else {
	print "\t\tCould not locate executable for $prog inside the package\n";
    }
    return $final;
}

sub test_perl {
    my ($prog, $path, $args, $out, @trash) = @_;
    my $final;
    if (-e $path) {
	print "\t\t$prog source found inside the package\n";
	# Use the absolute path of the executable
	my $cmd = "perl " . $ENV{"PWD"} . "/" . $path;

	# Check if it is working correctly
	$final = &test_cmd($cmd, $args, $out, "\t\t$prog (found in the package)", @trash);
    } else {
	print "\t\tCould not locate source for $prog inside the package\n";
    }
    return $final;

}

sub test_cmd {
    # Command, arguments, output, text, to be removed
    my ($cmd, $args, $out, $print, @trash) = @_;
    my $final;
    # Check if it is working correctly
    system("$cmd $args");
    if ($? == 0 && -s $out) {
	print "$print is working correctly\n";
	# Save the command that works
	$final = $cmd;
    } else {
	print "$print is not working correctly\n";
    }
    # Remove possible output files
    for (@trash) {
	unlink $_ if -e $_;
    }
    return $final;
}
