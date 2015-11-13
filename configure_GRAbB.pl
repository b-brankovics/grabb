#!/usr/bin/perl -w
use strict;
use autodie;


print "This is a script to test the installation of prerequisites for the programs in\n" .
    "the GRAbB package, and to configure the scripts before usage.\n";


print "Configuring GRAbB package\n";

# Variables that are true if given element is functional
my ($bait, $collect, $edena, $velveth, $velvetg, $exonerate);

# Read commands from GRAbB.pl
open(my $source, '<', 'GRAbB.pl') || die "Could not open GRAbB.pl. Make sure GRAbB.pl is in the base directory of the package and run configure_GRAbB.pl there.\n";
for(<$source>) {
    if (/my\s+\$bait_cmd\s+=\s+"([^"]+)"/) {
	$bait = $1;
    } elsif (/my\s+\$collect_cmd\s+=\s+"([^"]+)"/) {
	$collect = $1;
    } elsif (/my\s+\$edena_cmd\s+=\s+"([^"]+)"/) {
	$edena = $1;
    } elsif (/my\s+\$velveth_cmd\s+=\s+"([^"]+)"/) {
	$velveth = $1;
    } elsif (/my\s+\$velvetg_cmd\s+=\s+"([^"]+)"/) {
	$velvetg = $1;
    } elsif (/my\s+\$exonerate_cmd\s+=\s+"([^"]+)"/) {
	$exonerate = $1;
    }
}
close $source;

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

    $working = &test_cmd($bait, $call, $out, "\tRead baiting command found in GRAbB.pl source code", @trash) if $bait;

    $working = &test_path($prog, $call, $out, @trash) unless $working;
    $working = &test_exe($prog, $path, $call, $out, @trash) unless $working;

    print "\tConfiguration script failed to locate a working version of $prog\n" unless $working;
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

    $working = &test_cmd($collect, $call, $out, "\tRead collecting command found in GRAbB.pl source code", @trash) if $collect;

    $working = &test_path($prog, $call, $out, @trash) unless $working;

    print "\tConfiguration script failed to locate a working version of $prog\n" unless $working;
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
    if ($edena) {
	my $cmd = $edena;
	my $print = "\t\t$prog (found in the GRAbB.pl source)";
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
    }
    unless ($working) {
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
    }

    unless ($working) {    
	my $print = "\t\t$prog (found in the package)";
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
    if ($velveth) {
	# Test if the program is available through the path
	my $print = "\t\t$prog (found in the GRAbB.pl source)";
	my $cmd = $velveth;
	    
	# Check if it is working correctly
	system("$cmd $call");
	if ($? == 0 && -s $out) {
	    print "$print is working correctly\n";
	    # Save the command that works
	    $working = $cmd;
	} else {
	    print "$print is not working correctly\n";
	}
    }
    unless ($working) {
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
    }
    unless ($working) {    
	
	my $print = "\t\t$prog (found in the package)";
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
    if ($velvetg) {
	my $print = "\t\t$prog (found in the GRAbB.pl source)";
	my $cmd = $velvetg;
	    
	# Check if it is working correctly
	system("$cmd $call");
	if ($? == 0 && -s $out) {
	    print "$print is working correctly\n";
	    # Save the command that works
	    $working = $cmd;
	} else {
	    print "$print is not working correctly\n";
	}
    }
    unless ($working) {
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
    }
    unless ($working) {    
	
	my $print = "\t\t$prog (found in the package)";
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

    $working = &test_cmd($exonerate, $call, $out, "\tExonerate command found in GRAbB.pl source code", @trash) if $exonerate;

    $working = &test_path($prog, $call, $out, @trash) unless $working;

    print "\tConfiguration script failed to locate a working version of $prog\n" unless $working;

    $exonerate = $working;
}

# Create configured copy of GRAbB.pl in bin directory
my $summary = "\n" . "=" x 80 . "\n";
if ($bait && $collect) {
    my %pair;
    $pair{'(my\s+\$bait_cmd\s+=\s+")[^"]+"'} = $bait;
    $pair{'(my\s+\$collect_cmd\s+=\s+")[^"]+"'} = $collect;
    $pair{'(my\s+\$edena_cmd\s+=\s+")[^"]+"'} = $edena if $edena;
    $pair{'(my\s+\$velveth_cmd\s+=\s+")[^"]+"'} = $velveth if $velveth;
    $pair{'(my\s+\$velvetg_cmd\s+=\s+")[^"]+"'} = $velvetg if $velvetg;
    $pair{'(my\s+\$exonerate_cmd\s+=\s+")[^"]+"'} = $exonerate if $exonerate;
    &setup_program('GRAbB.pl', 'GRAbB.pl', \%pair);
    
    $summary .= "\nGRAbB.pl configuration:\n";
    
    if ($edena && $velveth && $velvetg && $exonerate) {
	$summary .= "\tGRAbB.pl is fully configured\n";
    } else {
	$summary .= "\tGRAbB.pl is minimally configured\n";
    }
    if ($edena) {
	$summary .= "\t\tEdena assembler is setup correctly\n";
    }
    if ($velveth && $velvetg) {
	$summary .= "\t\tVelvet assembler is setup correctly\n";
    }
    unless ($edena || ($velveth && $velvetg)) {
	$summary .= "\t\tGRAbB.pl cannot run without an external assembler '--assembler external <extranal_skeleton.pl>'\n"
    }
    if ($exonerate) {
	$summary .= "\t\tExonerate is setup correctly => exonerate mode can be used\n";
    } else {
	$summary .= "\t\tExonerate is not setup correctly => exonerate mode is disabled\n";
    }
} else {
    $summary .= "\nCould not configure GRAbB.pl\n";
    $summary .= "\tConfiguration script failed to locate a working read baiting program\n" unless $bait;
    $summary .= "\tConfiguration script failed to locate a working read collecting program\n" unless $collect;
}


# For PRINSEQ-lite the following modules are needed
{# Check helper scripts
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
	print "\tAll the Perl modules are installed that are needed to run PRINSEQ-lite.\n";
	# Configure helper scripts
	$summary .= "\nhelper_program configuration\n";
	&configure_helper('helper_programs/single2pairs', 'single2pairs');
	&configure_helper('helper_programs/uniform_length', 'uniform_length');
    }
}

print $summary;

# Copy rest of the helper programs
unless (-d "bin") {
    mkdir "bin" || die "Could not create bin directory\n";
    print "bin directory is created to store the configured executables\n"
}
my @rest = qw/fasta_shift  fastq2fasta  get_overlaps  interleaved2pairs  merge_contigs  pairwise_alignment_score  rename_fastq  reverse_complement/;
for (@rest) {
    my $old = "helper_programs/$_";
    my $new = "bin/$_";
    system("cp $old $new");
    chmod 0755, $new;
}

#==========Subroutines===============================================================
sub configure_helper {
    # Configure helper scripts
    my ($file, $name) = @_;
    my $prinseq;
    open(my $helper, '<', $file) || die "Could not open $file. Make sure yoe run configure_GRAbB.pl in the base directory of the package.\n";
    for (<$helper>) {
	if (/my\s+\$prinseq_cmd\s+=\s+"([^"]+)"/) {
	    $prinseq = $1;
	}
    }
    close $helper;

    my $working;
   
    if ($prinseq) {
	my $cmd = $prinseq;
	my $print = "\tprinseq-lite.pl command found in $file source";
	system("$cmd >/dev/null 2>/dev/null");
	if ($? == 0) {
	    print "$print is working correctly\n";
	    # Save the command that works
	    $working = $cmd;
	} else {
	    print "$print is not working correctly\n";
	}
    }

    unless ($working) {
	my $cmd = "prinseq-lite.pl";
	my $print = "\t\tprinseq-lite.pl found in the path";
	if (`which $cmd`) {
	    print "\t\tprinseq-lite found in the path\n";
	    system("$cmd >/dev/null 2>/dev/null");
	    if ($? == 0) {
		print "$print is working correctly\n";
		$working = $cmd;
	    } else {
		print "$print is not working correctly\n";
	    }
	} else {
	    print "$cmd is not in the path\n";
	    my $print = "prinseq-lite executable";
	    $cmd = $ENV{"PWD"} . '/' . '3rd_party_programs/prinseq-lite-0.20.4/prinseq-lite.pl';
	    if (-x $cmd) {
		print "\t\tprinseq-lite executable is found in the package\n";
		system("$cmd >/dev/null 2>/dev/null");
		if ($? == 0) {
		    print "$print is working correctly\n";
		    $working = $cmd;
		} else {
		    print "$print is not working correctly\n";
		}
	    } else {
		print "\t\tprinseq-lite executable is not found in the package\n";
		if (-s $cmd) {
		    $cmd = "perl $cmd";
		    print "\t\tprinseq-lite source is found in the package\n";
		    $print = "prinseq-lite source";
		    system("$cmd >/dev/null 2>/dev/null");
		    if ($? == 0) {
			print "$print is working correctly\n";
			$working = $cmd;
		    } else {
			print "$print is not working correctly\n";
		    }
		} else {
		    print "\t\tprinseq-lite source is not found in the package\n";
		}
	    }
	}
    }
    $prinseq = $working;
    if ($prinseq) {
	my %pair;
	$pair{'(my\s+\$prinseq_cmd\s+=\s+")[^"]+"'} = $prinseq;
	&setup_program($file, $name, \%pair)
    }
}

sub setup_program {
    my ($read, $name, $replace) = @_;
    # Put the final program to bin directory inside GRAbB package
    my $new = "bin/$name";
    unless (-d "bin") {
	mkdir "bin" || die "Could not create bin directory\n";
	print "bin directory is created to store the configured executables\n"
    }
    open(my $in, '<', $read) || die $!;
    open(my $out, '>', $new) || die $!;
    for (<$in>) {
	for my $pattern (keys %{$replace}) {
	    if (/$pattern/) {
		s/$pattern/$1$replace->{$pattern}\"/;
	    }
	}
	print {$out} $_;
    }
    close $in;
    close $out;
    print "Configured $name can be found in the bin folder ($new)\n";
    chmod 0755, $new || die $!;
    print "Configured $name can be found in the bin folder ($new) is made executable\n";
    $summary .= "\t$name is configured\n" unless $name eq 'GRAbB.pl';
}

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
