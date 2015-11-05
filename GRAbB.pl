#!/usr/bin/perl
use strict;
use File::Spec;
use POSIX qw(strftime);
use autodie;

#======================================HEADER================================================================================================

##############################################################################################################
##############################################################################################################
##                                                                                                          ##
##############################################################################################################
##############################################################################################################
## This program let's you selcetively assemble parts of a genome from raw reads. The program works in an    ##
## iterative fashon: finds reads that match to bait/reference using exact kmer match, assembles the reads   ##
## de novo.                                                                                                 ##
## For more information read the documentation or the usage message bellow (also by invoking the program    ##
## without any arguments.                                                                                   ##
##############################################################################################################
##############################################################################################################

#======================================USAGE=================================================================================================

#######################################################################################################################
#                                                        USAGE                                                        #
#######################################################################################################################
# Usage message to be printed for incorect invocation                                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
my $usage = "Usage:\n" .                                                                                              #
"$0 --ref <reference file> --reads <read file 1> [<read file 2>] --folder <directory> --prefix <prefix> [options]\n" .#
"\n" .                                                                                                                #
"\t--ref <ref.fas>\n" .                                                                                               #
"\t\t\t\tReference file: FASTA format, the id lines must be unique\n" .                                               #
"\t\t\t\t(only the first word is used as an id)\n" .                                                                  #
"\t--reads	<r1.fastq> [<r2.fastq> ...]\n" .                                                                      #
"\t\t\t\tRead files in fasta or fastq format, may be compressed using\n" .                                            #
"\t\t\t\tgzip. If the reads are paired-end specify separately both the\n" .                                           #
"\t\t\t\tforward and the reverse read files. Otherwise specify only one\n" .                                          #
"\t\t\t\tfile or use '--single' option\n" .                                                                           #
"\t--folder <folder_name>\n" .                                                                                        #
"\t\t\t\tThe directory where all the output will be saved. If the\n" .                                                #
"\t\t\t\tdirectory is non-empty then the files it contains can be used\n" .                                           #
"\t\t\t\tlike internal files. In this manner previous runs can be continued,\n" .                                     #
"\t\t\t\tmake sure to remove or replace files that would suggest completion\n" .                                      #
"\t--prefix <prefix_of_output>\n" .                                                                                   #
"\t\t\t\tThe prefix for the output files\n" .                                                                         #
"Additional options:\n" .                                                                                             #
"\t--single\n" .                                                                                                      #
"\t\t\t\tTreat reads as unpaired reads even if two read files are specified\n" .                                      #
"\t--bait <bait.fas>\n" .                                                                                             #
"\t\t\t\tA separate bait file can be specified besides the reference file,\n" .                                       #
"\t\t\t\tthis file together with the reference file will be used as first bait.\n" .                                  #
"\t\t\t\tUseful when using special criterion for the assembly, such as homology\n" .                                  #
"\t--min_length=<int>\n" .                                                                                            #
"\t\t\t\tMinimum size required for a contig to be included for completion testing and baiting\n" .                    #
"\t--type <run_type_string>\n" .                                                                                      #
"\t\t\t\tSpecial critera can be specified for the run\n" .                                                            #
"\t\t\t\t\tWhen should the iteration stop?\n" .                                                                       #
"\t\t\t\t\t\tThe intrinsic criteria are:\n" .                                                                         #
"\t\t\t\t\t\t\tThere are no new reads\n" .                                                                            #
"\t\t\t\t\t\t\tThere is no assembly\n" .                                                                              #
"\t\t\t\t\t\t\tThe bait sequence did not change\n" .                                                                  #
"\t\t\t\t\t\tExtra criteria:\n" .                                                                                     #
"\t\t\t\t\t\t\tLength of the assembly\n" .                                                                            #
"\t\t\t\t\t\t\tLength of the longest contig of the assembly\n" .                                                      #
"\t\t\t\t\t\t\tN50 value of the assembly\n" .                                                                         #
"\t\t\t\t\t\t\tReference matches to assembly (uses exonerate)\n" .                                                    #
"\t\t\t\t\tOr to have parallel runs (threads)\n" .                                                                    #
"\t\t\t\t\t\tPossible values:\n" .                                                                                    #
"\t\t\t\t\t\tmulti		Activate multi mode (parallel runs)\n" .                                              #
"\t\t\t\t\t\ttotal=<int>	Minimum size (bp) of the assembly\n" .                                                #
"\t\t\t\t\t\tlongest=<int>	Minimum size (bp) of the longest contig\n" .                                          #
"\t\t\t\t\t\tN50=<int>	Minimum N50 value (bp) of the assembly\n" .                                                   #
"\t\t\t\t\t\texonerate	Matches reference using exonerate.\n" .                                                       #
"\t\t\t\t\t\t\t\t(Matched region is saved separately)\n" .                                                            #
"\t\t\t\t\t\tdiff            Uses the identifier lines from the reference\n" .                                        #
"\t\t\t\t\t\t                to get the completion criterion for each thread\n" .                                     #
"\t\t\t\t\t\t                (These criteria overwrite the globally defined\n" .                                      #
"\t\t\t\t\t\t                criteria. Write \"exhaustive\" in the identifier\n" .                                    #
"\t\t\t\t\t\t                line if exhaustive should be used.)\n" .                                                 #
"\t\t\t\t\t\t                (Need to be used together with multi)\n" .                                               #
"\t\t\t\tMultiple options can be used at the same time, but then they have to be\n" .                                 #
"\t\t\t\ttyped as a concatenated string (or using quotation)\n" .                                                     #
"\t\t\t\t(e.g. --type <multiexonerate | multi_exonerate | \"multi exonerate\">)\n" .                                  #
"\t--arg1 \"argument1 argument2 argument3\"\n" .                                                                      #
"\t\t\t\tArguments needed for the assembler program used for graph/hash generation\n" .                               #
"\t--arg2 \"argument1 argument2 argument3\"\n" .                                                                      #
"\t\t\t\tArguments needed for the assembler program used for the assembly\n" .                                        #
"\t--assembler assembler_name\n" .                                                                                    #
"\t\t\t\tSpecify the assembler to be used: edena, velvet, alternative or external\n" .                                #
"\t--clean\n" .                                                                                                       #
"\t\t\t\tRemove some internal files to save disk space.\n" .                                                          #
"\t\t\t\tIf specified twice, then only result files are kept, the rest is deleted\n";                                 #
#######################################################################################################################

#======================================START=UP==============================================================================================

####################################################################################################
####################################################################################################
##                                     VARIABLES                                                  ##
####################################################################################################
####################################################################################################
# A scalar variable to store the log information before log file is created                        #
my $p_log;                                                                                         #
# A scalar variable to store the print message before printing it to STDOUT and the log file       #
#	Didn't use IO:Tee, because it is not part of the basic Perl installation)                  #
my $print;                                                                                         #
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
# Constant variables                                                                               #
#	Only declared at the start and never changed afterwards                                    #
#	These variables have to be adjusted to the environment                                     #
my $bait_cmd =    "mirabait";       # The command to invoke the baiting program                    #	!!!!
my $collect_cmd = "seqtk subseq";   # The command to invoke the read collecting program            #	!!!!
my $edena_cmd =	  "edena";          # The command to invoke the default assembler program          #	!!!!
my $velveth_cmd = "velveth";        # The command to invoke velvet graph/hash generator program    #	!!!!
my $velvetg_cmd = "velvetg";        # The command to invoke velvet assembler program               #	!!!!
my $alt_1_cmd =   "";               # The command to invoke graph/hash generator program           #	!!!!!!	alternative assembler (graph)
my $alt_2_cmd =   "";               # The command to invoke assembler program                      #	!!!!!!	alternative assembler (assembly)
my $ext_assembler = "";             # A perl script that does the assembly                         #
my $exonerate_cmd = "exonerate";    # The command to invoke exonerate                              #	!!!!
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
# The variables holding file names and / or depend on the command line invocation                  #
my $ref;            # What is the reference file?                                                  #
my @reads;          # What are the reads files used as input?                                      #
my $format;         # What is the file format of the read files? (fasta or fastq)                  #
my $gzip;           # True if the reads are zipped using gzip                                      #
my $single;         # If reads are paired FALSE else TRUE                                          #
my $folder;         # What is the working directory?                                               #
my $prefix;         # What is the prefix for the output?                                           #
my $type = "";      # What is the criteria to stop?  If undef then when there are no new reads.    #
my $min_len;        # What is the minimal length for contigs to be used?                           #
my @parameters;     # What are the parameters for the assembler?                                   #
my $extra_bait;     # Was there extra bait specified when launching in multi mode?                 #
my $clean = 0;      # 1 if some files should be deleted, 2 if only results should be kept          #
my $log;            # The file handle for the log file                                             #
my $bait;           # What is the bait file?                                                       #
my $old_collection; # The collection file of the previous iteration (first iteration = null file)  #
my $new_collection; # The collection file of the current iteration                                 #
my @readpools;      # The read files generated by baiting and used by assembly                     #
my $assembly_file;  # The FASTA file containing the assembly                                       #
my @current_asmbls; # The list of assembly files to be modified if min_len is set                  #
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
# Parse the commandline input and load the variables                                               #
#   This also creates the folder and all the files that are based on the input                     #
#       (reference, bait, links to the read files)                                                 #
&parse_input(\@ARGV, \$ref, \@reads, \$single, \$folder, \$prefix, \$type, \@parameters, \$p_log,  #
             \$bait, \$extra_bait, \$clean, \$gzip, \$min_len, \$ext_assembler, \$format);         #
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
# Test if all the specified programs are usable                                                    #
for ($bait_cmd, $collect_cmd) {                                                                    #
    die "$_ program is not executable\n" unless `which $_`;                                        #
}                                                                                                  #
# Test if the assembler is usable or die                                                           #
if ($parameters[0] eq "alternative") {                                                             #
    die "$alt_1_cmd program is not executable\n" unless `which $alt_1_cmd`;                        #
    die "$alt_2_cmd program is not executable\n" unless `which $alt_2_cmd`;                        #
} elsif ($parameters[0] eq "velvet") {                                                             #
    die "$velveth_cmd program is not executable\n" unless `which $velveth_cmd`;                    #
    die "$velvetg_cmd program is not executable\n" unless `which $velvetg_cmd`;                    #
} else {                                                                                           #
    die "$edena_cmd program is not executable\n" unless `which $edena_cmd`;                        #
}                                                                                                  #
# If exonerate is necessary then test it also                                                      #
if ($type =~ /exonerate/) {                                                                        #
    die "$exonerate_cmd program is not executable\n" unless `which $exonerate_cmd`;                #
}                                                                                                  #
####################################################################################################
# Dynamic variables (They are undefined at start)                                                  #
my $done;           # Is the test criteria met?                                                    #
my $round;          # The iteration counter                                                        #
my $current_folder; # The name of the current round folder                                         #
my $previous_folder;# The name of the previous round folder                                        #
my %multi;          # Contains all the threads: keys: thread id; values: status (undef if active)  #
#                   #   There are subroutines that operate on this hash                            #
#                   #       load_multi      creates the threads at the start                       #
#                   #       get_actives     returns all the threads in a sorted order              #
my %types;          # Contains all the threads: keys: thread id; values: criterion type            #
####################################################################################################
# Declare some more variables                                                                      #
#   Load the file names that do not depend on the command line input                               #
$assembly_file = "assembly.fas";                                                                   #
$new_collection = "new_collection.list";                                                           #
$old_collection = "old_collection.list";                                                           #
{   # Using a bare block to create local variables that will not interfere later                   #
    # Create old_collection file if it does not exist already                                      #
    #   Using append to open the file so content is not destroyed                                  #
    open my $temp_handle, '>>', "$folder\/old_collection.list";                                    #
    close $temp_handle;                                                                            #
    # For each read file there should be a readpool file used                                      #
    # Create a counter                                                                             #
    my $i;                                                                                         #
    for (@reads) {                                                                                 #
        $i++;                                                                                      #
        # Save the relative path of the readpool file into the readpool array                      #
        push @readpools, "readpool$i\.$format";                                                    #
    }                                                                                              #
}                                                                                                  #
####################################################################################################

#=========================================MAIN=PART==========================================================================================

###############################################################################################
#                                         PREPARATION                                         #
###############################################################################################			Working directory
# Move to working directory and create more variables and files                               #			(the specified folder is marked "/")
chdir ($folder);                                                                              #			/
# Open log file for appending, even if it already exists                                      #			/
open $log, '>>', "$prefix\.log";                                                              #			/
                                                                                              #			/
# Print greeting text                                                                         #			/
$print = "\n" . "X" x 90 . "\n";                                                              #			/
print {$log} $print; print $print;                                                            #			/
$print = "GRAbB: Genome Region Assembly by Baiting\nInvocation:\t$0 @ARGV\n";                 #			/
print {$log} $print; print $print;                                                            #			/
# Print what belongs into the log                                                             #			/
print {$log} $p_log; print $p_log;                                                            #			/
                                                                                              #			/
&load_multi(\$log, \$ref, \$bait, \$type, \%multi, \%types);                                  #			/
                                                                                              #			/
# Detect previous runs to be continued or start a new one                                     #			/
{                                                                                             #			/
    my @list = glob("Round*");                                                                #			/
                                                                                              #			/
    # Find the "Round" folder with the highest number and continue the run from that cycle    #			/
    my ($highest) = sort {$b <=> $a} map {/^Round(\d+)/; $1} @list;                           #			/
                                                                                              #			/
    # First step of the loop is incrementing the round variable                               #			/
    if ($highest) {                                                                           #			/
        $print = "Detected previous run, the folders and files will be used\n";               #			/
        print {$log} $print; print $print;                                                    #			/
        $round = $highest - 1;                                                                #			/
    } else {                                                                                  #			/
        $round = 0;                                                                           #			/
    }                                                                                         #			/
}                                                                                             #			/
###############################################################################################			/
####################################################################################################	/
# Main loop of the program                                                                         #	/
while (not $done) {                                                                                #	/
    $round++; # increment the counter                                                              #	/
    # Save the name of the folder name of the previous round                                       #	/
    $previous_folder = $current_folder;                                                            #	/
    # Compile the name of the folder of the current round                                          #	/
    if (length $round > 3) {                                                                       #	/
        $current_folder = "Round$round";                                                           #	/
    } else {                                                                                       #	/
        $current_folder = "Round" . ("0" x (3 - length $round)) . $round;                          #	/
    }                                                                                              #	/
                                                                                                   #	/
    # Print a round separator line                                                                 #	/
    $print = "-" x 90 . "\n";                                                                      #	/
    print {$log} $print; print $print;                                                             #	/
                                                                                                   #	/
    # Print time and the current round                                                             #	/
    $current_folder =~ /Round(\d+)/;                                                               #	/
    $print = "" . strftime('%H:%M:%S',localtime) . "\tRound $1\n";                                 #	/
    print {$log} $print; print $print;                                                             #	/
                                                                                                   #	/
    # Create the current round folder if it is not present                                         #	/
    unless (-d $current_folder) {                                                                  #	/
        mkdir $current_folder;                                                                     #	/
    }                                                                                              #	/
                                                                                                   #	/
    # Move to the current round folder                                                             #	/
    chdir $current_folder;                                                                         #	/					=>	/Round###
    # Update the relative paths of the files that are in the main directory                        #	/Round###
    #   The file names are fixed inside the subroutine                                             #	/Round###
    #   The prefix is passed then the reference file, bait file and old collection variables       #	/Round###
    &update_path(\"../", \$ref, \$bait, \$old_collection);#/                                       #	/Round###
                                                                                                   #	/Round###
    # Check if all the necessary files are in order                                                #	/Round###
    #   Creates a folder called test_reads which will be destroyed if all goes well                #	/Round###
    #   The program dies if some error is found                                                    #	/Round###
    &check_files(\$log, \$bait, \$ref, \@reads, \$single, \$gzip);                                 #	/Round###
                                                                                                   #	/Round###
    # Find reads and create read(pool) file(s) for the assembly                                    #	/Round###
    #   Stop if there are no new reads (Returns "No new reads")                                    #	/Round###
    $done = &find_reads(\$log, \$bait, \@reads, \$old_collection, \$new_collection, \@readpools);  #	/Round###
    # If done then break the loop and proceed to the finishing steps                               #	/Round###
    if ($done) {                                                                                   #	/Round###
        chdir "..";                                                                                #	/Round###			=>	/
        last;                                                                                      #==================== Breaks the Main loop
    }                                                                                              #
                                                                                                   #	/Round###
    # If clean mode is selected remove the readpool of the previous round                          #	/Round###
    #   If there was a previous round already                                                      #	/Round###
    #   The new readpool files contain all the reads that were in the previous ones                #	/Round###
    if ($clean && $previous_folder) {                                                              #	/Round###
        for (@readpools) {                                                                         #	/Round###
            system("rm", "../$previous_folder/$_") if -e "../$previous_folder/$_";                 #	/Round###
        }                                                                                          #	/Round###
    }                                                                                              #	/Round###
                                                                                                   #	/Round###
    # Copy the new collection to the old collection                                                #	/Round###
    system("cp", $new_collection, $old_collection);                                                #	/Round###
                                                                                                   #	/Round###
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ##	/Round###
    # Run a separate analysis for each of the threads                                              #	/Round###
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ##	/Round###
    # If multi mode is selected then does a separate run for each thread of the active threads     #	/Round###
    # Otherwise, only runs on the single main thread                                               #	/Round###
    for my $thread (&get_actives(\%multi)) {                                                       #	/Round###
        # Create the current thread folder if it is not present                                    #	/Round###
        unless (-d $thread) {                                                                      #	/Round###
            mkdir $thread;                                                                         #	/Round###
        }                                                                                          #	/Round###
                                                                                                   #	/Round###
        # Move to the current thread folder                                                        #	/Round###
        chdir $thread;                                                                             #	/Round###			=>	/Round###/thread#
                                                                                                   #	/Round###/thread#
        # Create the local "links" to the necessary files                                          #	/Round###/thread#
        #   Use a separate name for the thread specific files                                      #	/Round###/thread#
        my ($thread_ref, $thread_bait, $thread_old) = ($ref, $bait, $old_collection);              #	/Round###/thread#
                                                                                                   #	/Round###/thread#
        # Update the relative paths of the files that are in the main thread directory             #	/Round###/thread#
        #   The file names are fixed inside the subroutine                                         #	/Round###/thread#
        #   The prefix is passed then the reference file, bait file and old collection variables   #	/Round###/thread#
        &update_path(\"../../$thread\/", \$thread_ref, \$thread_bait, \$thread_old);#"             #	/Round###/thread#
                                                                                                   #	/Round###/thread#
        # Create the symlinks to the readpool files of the Round folder so that the same read file #	/Round###/thread#
        #   name can be used for the threads                                                       #	/Round###/thread#
        my $i = 0;                                                                                 #	/Round###/thread#
        for (@readpools) {                                                                         #	/Round###/thread#
            $i++;                                                                                  #	/Round###/thread#
            unless (-l "../reads$i\.$format") {                                                    #	/Round###/thread#
                symlink(File::Spec->rel2abs("../$_"), "../reads$i\.$format");                      #	/Round###/thread#
            }                                                                                      #	/Round###/thread#
        }                                                                                          #	/Round###/thread#
                                                                                                   #	/Round###/thread#
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ##	/Round###/thread#
        # An inner loop for the individual thread in multi mode with extra bait at startup         #	/Round###/thread#
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ##	/Round###/thread#
        # A variable that shows if another iteration step is needed or not                         #	/Round###/thread#
        my $continue = 1;                                                                          #	/Round###/thread#
        # A counter for the number of the iteration                                                #	/Round###/thread#
        my $iteration;                                                                             #	/Round###/thread#
        while ($continue > 0) {                                                                    #	/Round###/thread#
            # Default is to only have one cycle => basically not a loop                            #	/Round###/thread#
            $continue--;                                                                           #	/Round###/thread#
            $iteration++;                                                                          #	/Round###/thread#
                                                                                                   #	/Round###/thread#
            # Check if all the necessary files are in order                                        #	/Round###/thread#
            &check_files(\$log, \$thread_bait, \$thread_ref, \@reads, \$single, \undef);           #	/Round###/thread#
                                                                                                   #	/Round###/thread#
            # Find reads and create read(pool) file(s) for the assembly                            #	/Round###/thread#
            #   Stop if there are no new reads (Returns "No new reads")                            #	/Round###/thread#
            #   Only needed in multi mode                                                          #	/Round###/thread#
            if ($type && $type =~ /multi/) {                                                       #	/Round###/thread#
                # Print thread separator line                                                      #	/Round###/thread#
                $print = "\t\t\t$thread\t" . ("-" x 50) . "\n";                                    #	/Round###/thread#
                print {$log} $print; print $print;                                                 #	/Round###/thread#
                                                                                                   #	/Round###/thread#
                # Find reads for the current thread and update the threads state                   #	/Round###/thread#
                $multi{$thread} = &find_reads(\$log, \$thread_bait, \@reads,                       #	/Round###/thread#
                                              \$thread_old, \$new_collection,                      #	/Round###/thread#
                                              \@readpools, \"\t"); #"                              #	/Round###/thread#
                                                                                                   #	/Round###/thread#
                # If no new reads were found then go to the next thread                            #	/Round###/thread#
                if ($multi{$thread}) {                                                             #	/Round###/thread#
                    # If this is the first Round with multi mode and extra bait,                   #	/Round###/thread#
                    #	Then mark it as active thread                                              #	/Round###/thread#
                    if ($extra_bait) {                                                             #	/Round###/thread#
                        # During the previous inner iteration (if there was one) the files created #	/Round###/thread#
                        #   were renamed by adding the "old_" prefix. Because there was no new     #	/Round###/thread#
                        #   information during this inner iteration the old files are renamed ( by #	/Round###/thread#
                        #   removing the "old_" prefix) and overwrite the current iteration files. #	/Round###/thread#
                        for (glob("old_*")) {                                                      #	/Round###/thread#
                            my $rest = $_;                                                         #	/Round###/thread#
                            $rest =~ s/old_//;                                                     #	/Round###/thread#
                            system("mv", $_, $rest);                                               #	/Round###/thread#
                        }                                                                          #	/Round###/thread#
                        # If there were no reads found then stop the thread                        #	/Round###/thread#
                        $multi{$thread} = undef unless $multi{$thread} eq "No reads";              #	/Round###/thread#
                    }                                                                              #	/Round###/thread#
                    # Because $continue is zero, there won't be a next inner loop                  #	/Round###/thread#
                    next;                                                                          #==================== Breaks the Inner loop
                }                                                                                  #
            } else {                                                                               #	/Round###/thread#
                # IF NOT in multi mode                                                             #	/Round###/thread#
                #   THEN create symlinks to the selected reads in the thread folder                #	/Round###/thread#
                #       Because the find reads subroutine was skipped                              #	/Round###/thread#
                for (@readpools) {                                                                 #	/Round###/thread#
                    symlink(File::Spec->rel2abs("../$_"), "$_") unless -e $_;                      #	/Round###/thread#
                }                                                                                  #	/Round###/thread#
                # And copy the new collection file into the current working directory              #	/Round###/thread#
                system("cp", "../$new_collection", $new_collection);                               #	/Round###/thread#
            }                                                                                      #	/Round###/thread#
                                                                                                   #	/Round###/thread#
            # Do assembly according to parameters and create the assembly file                     #	/Round###/thread#
            &assemble(\$log, \@readpools, \$single, \@parameters,                                  #	/Round###/thread#
                      \$assembly_file, \"\t", \$format);#"                                         #	/Round###/thread#
                                                                                                   #	/Round###/thread#
            # Copy the latest assembly to the thread folder in the main directory                  #	/Round###/thread#
            system("cp", $assembly_file, "../../$thread/$assembly_file");                          #	/Round###/thread#
                                                                                                   #	/Round###/thread#
            if ($clean < 2) {                                                                      #	/Round###/thread#
                # Unless hyperclean mode is selected, keep a copy of each assembly file            #	/Round###/thread#
                #   Get the next unused serial number for the assembly file                        #	/Round###/thread#
                #       The first serial number should be 1                                        #	/Round###/thread#
                my $next = 1;                                                                      #	/Round###/thread#
                while (-e "../../$thread/assembly_$next\.fas") {                                   #	/Round###/thread#
                    $next++;                                                                       #	/Round###/thread#
                }                                                                                  #	/Round###/thread#
                # copy the assembly to the new file                                                #	/Round###/thread#
                system("cp", "../../$thread/$assembly_file", "../../$thread/assembly_$next\.fas"); #	/Round###/thread#
                if ($min_len) {                                                                    #	/Round###/thread#
                    system("cp", "../../$thread/$assembly_file",                                   #	/Round###/thread#
                           "../../$thread/assembly_$next\_all\.fas");                              #	/Round###/thread#
                    @current_asmbls = ("../../$thread/$assembly_file",                             #	/Round###/thread#
                                       "../../$thread/assembly_$next\.fas");                       #	/Round###/thread#
                }                                                                                  #	/Round###/thread#
            }                                                                                      #	/Round###/thread#
                                                                                                   #	/Round###/thread#
            # Check if another iteration is needed or not                                          #	/Round###/thread#
            #   There are different types of criteria                                              #	/Round###/thread#
            #       constant                                                                       #	/Round###/thread#
            #           no new reads           =>   stop (already checked before the assembly)     #	/Round###/thread#
            #           bait hasn't changed    =>	stop (returns "Same bait")                 #	/Round###/thread#
            #                                   (in the next iteration there would be no new reads)#	/Round###/thread#
            #           assembly file is empty =>   stop (returns "No assembly")                   #	/Round###/thread#
            #       optional                                                                       #	/Round###/thread#
            #           length of the assembly                                                     #	/Round###/thread#
            #           "completeness" using homology                                              #	/Round###/thread#
            #   If not done then replace the old bait file with the new bait                       #	/Round###/thread#
            if ($type && $type =~ /diff/) {                                                        #	/Round###/thread#
                $multi{$thread} = &test_completion(\$log, \$thread_ref, \$thread_bait,             #	/Round###/thread#
                                                    \$types{$thread}, \$assembly_file, \$min_len,  #	/Round###/thread#
                                                    \"\t", \@current_asmbls);#"                    #	/Round###/thread#
            } else {                                                                               #	/Round###/thread#
                $multi{$thread} = &test_completion(\$log, \$thread_ref, \$thread_bait, \$type,     #	/Round###/thread#
                                                   \$assembly_file, \$min_len, \"\t",#"            #	/Round###/thread#
                                                   \@current_asmbls);                              #	/Round###/thread#
            }                                                                                      #	/Round###/thread#
            # If homology is used as test there might be a "result file created that holds the     #	/Round###/thread#
            #   the sequence matching the reference from the assembly                              #	/Round###/thread#
            #   Copy this result file to the thread folder                                         #	/Round###/thread#
            system("cp", "result.fas", "../../$thread/result.fas") if -e "result.fas";             #	/Round###/thread#
                                                                                                   #	/Round###/thread#
            # Copy the new collection to the old collection                                        #	/Round###/thread#
            system("cp", $new_collection, $thread_old);                                            #	/Round###/thread#
                                                                                                   #	/Round###/thread#
            # If this is the first round of a multi run with extra bait specified                  #	/Round###/thread#
            #   Then rename the files in the thread folder by adding the "old_" prefix             #	/Round###/thread#
            if ($extra_bait) {                                                                     #	/Round###/thread#
                # If there are files with "old_" prefix in the folder then delete them,            #	/Round###/thread#
                #   because all the files are present with updated data                            #	/Round###/thread#
                for (glob("old_*")) {                                                              #	/Round###/thread#
                    system("rm", $_);                                                              #	/Round###/thread#
                }                                                                                  #	/Round###/thread#
                                                                                                   #	/Round###/thread#
                # If the current thread is still active do a new inner iteration                   #	/Round###/thread#
                if (not $multi{$thread}){                                                          #	/Round###/thread#
                    # This means there will be a new iteration                                     #	/Round###/thread#
                    $continue++;                                                                   #	/Round###/thread#
                    # Rename current files by adding the "old_" prefix,                            #	/Round###/thread#
                    #	so that they can be restored later if no new reads are found               #	/Round###/thread#
                    for (glob("*")) {                                                              #	/Round###/thread#
                        system("mv", $_, "old_$_");                                                #	/Round###/thread#
                    }                                                                              #	/Round###/thread#
                } else {                                                                           #	/Round###/thread#
                    # If the thread terminated with "Same bait" value after more than one inner    #	/Round###/thread#
                    #   loop => Then reactivate the thread                                         #	/Round###/thread#
                    if ($multi{$thread} eq "Same bait" || $multi{$thread} eq "No improvement") {   #	/Round###/thread#
                        $multi{$thread} = undef if $iteration > 1;                                 #	/Round###/thread#
                    }                                                                              #	/Round###/thread#
                }                                                                                  #	/Round###/thread#
            }                                                                                      #	/Round###/thread#
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ##	/Round###/thread#
        # The end of the inner loop for the current thread                                         #	/Round###/thread#
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ##	/Round###/thread#
        }                                                                                          #	/Round###/thread#
                                                                                                   #	/Round###/thread#
        # Move back to the Round directory, because the first step of the thread loop is to move   #	/Round###/thread#
        #   into the thread directory.                                                             #	/Round###/thread#
        #   This way the working directory is the same as it was before the for and while loops    #	/Round###/thread#
        chdir "..";                                                                                #	/Round###/thread#	=>	/Round###/
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ##	/Round###
    # End of the (for) thread loop                                                                 #	/Round###
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ##	/Round###
    }                                                                                              #	/Round###
                                                                                                   #	/Round###
    # If clean mode is active, delete thread folders from Round folders which are not valuable     #	/Round###
    if ($clean && $previous_folder) {                                                              #	/Round###
        for my $thread (glob("thread_*")) {                                                        #	/Round###
            # Check which threads were active this turn                                            #	/Round###
            if ($multi{$thread} && $multi{$thread} eq "No new reads") {                            #	/Round###
                # If no new reads then there is nothing important in the current folder            #	/Round###
                system("rm", "-r", $thread);                                                       #	/Round###
            } else {                                                                               #	/Round###
                # The information content of the previous folder is a subset of the current's      #	/Round###
                system("rm", "-r", "../$previous_folder/$thread");                                 #	/Round###
            }                                                                                      #	/Round###
        }                                                                                          #	/Round###
    }                                                                                              #	/Round###
                                                                                                   #	/Round###
    # Move back to the parent (main) directory                                                     #	/Round###
    chdir "..";                                                                                    #	/Round###			=>	/
                                                                                                   #	/
    # Update the relative paths of the files that are in the main directory                        #	/
    #   The file names are fixed inside the subroutine                                             #	/
    #	The prefix is passed then the reference file, bait file and old collection variables       #	/
    &update_path(\"", \$ref, \$bait, \$old_collection);#"                                          #	/
                                                                                                   #	/
    # If this was the first round of a multi run with extra bait specified then make the           #	/
    #   $extra_bait variable undefined, because no more inner iterations will be need in the future#	/
    if ($extra_bait) {                                                                             #	/
        $extra_bait = undef;                                                                       #	/
    }                                                                                              #	/
                                                                                                   #	/
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ##	/
    # Prepare files for the next round if it is needed                                             #	/
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ##	/
    # Stop the main loop if there are no active threads left                                       #	/
    last if 0 == scalar( &get_actives(\%multi) );                                                  #==================== Breaks the Main loop
                                                                                                   #	/
    # Merge the individual bait files of the threads into a single bait file                       #	/
    #   Collect all the baits needed to be merged                                                  #	/
    my @baits = map{"$_/bait.fas"} &get_actives(\%multi);                                          #	/
    #   Merge the collected files and rename the entries so each one is uniq                       #	/
    `cat @baits | perl -ne 'if (/^>/) {\$n++; print ">\$n\n"} else {print}' >$bait`;               #	/
}                                                                                                  #	/
####################################################################################################

# Print separator line at the end of the main loop
$print = "=" x 90 . "\n";
print {$log} $print; print $print;

###############################################################################################################
# Final steps                                                                                                 #
###############################################################################################################
# Print the status of the threads and create the output files and do the necessary cleaning steps             # Final steps
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # Final steps
# Iterate through each thread                                                                                 # Final steps
for my $thread (sort{   my ($a1) = $a =~ /thread_(\d+)/;                                                      # Final steps
                        my ($b1) = $b =~ /thread_(\d+)/;                                                      # Final steps
                        $a1 <=> $b1                        } keys %multi) {                                   # Final steps
                                                                                                              # Final steps
    # If the thread was still active when the main loop stopped the reason of the stop should be copied       # Final steps
    #   to the thread status                                                                                  # Final steps
    $multi{$thread} = "$done" unless $multi{$thread};                                                         # Final steps
    # Possible values                                                                                         # Final steps
    #  Basic values                                                                                           # Final steps
    #   No assembly         => There was some error with the assembly (could be that there were to few reads) # Final steps
    #   Same bait           => The assembly was identical to the current bait                                 # Final steps
    #   No new reads        => There were no new reads found, thus the assembly would be the same             # Final steps
    #  Exonerate values                                                                                       # Final steps
    #   No improvement      => The exonerate analysis had the same result as during the previous round        # Final steps
    #   Multiple contigs    => The reference was covered by multiple contigs, the matched region won't be     # Final steps
    #                           saved to a separate file                                                      # Final steps
    #   Found               => The reference was completely matched to a single contig of the assembly        # Final steps
    #                           The matched region will be saved to a separate file                           # Final steps
    #  Length value                                                                                           # Final steps
    #   Reached:            => The specified minimal length was attained                                      # Final steps
                                                                                                              # Final steps
    # Print the threads final status                                                                          # Final steps
    $print = "$thread:\t$multi{$thread}\n";                                                                   # Final steps
    print {$log} $print; print $print;                                                                        # Final steps
                                                                                                              # Final steps
    # If there is a non-empty assembly file in the thread folder than copy it to the main folder              # Final steps
    if (-e "$thread/$assembly_file" && not -z "$thread/$assembly_file") {                                     # Final steps
        system("cp", "$thread/$assembly_file", "$prefix\_assembly_$thread\.fas");                             # Final steps
        # Also create a final assembly file in the thread folder                                              # Final steps
        system("cp", "$thread/$assembly_file", "$thread/final_assembly.fas");                                 # Final steps
                                                                                                              # Final steps
        # Print the location of the assembly output file                                                      # Final steps
        $print = "\tThe final assembly is $folder/$prefix\_assembly_$thread\.fas\n";                          # Final steps
        print {$log} $print; print $print;                                                                    # Final steps
    } else {                                                                                                  # Final steps
        # Print that there is no assembly file for the thread                                                 # Final steps
        $print = "\tThere is no assembly for this thread\n";                                                  # Final steps
        print {$log} $print; print $print;                                                                    # Final steps
    }                                                                                                         # Final steps
                                                                                                              # Final steps
    # If there is a result file in the thread folder then copy it to the main folder as an output file        # Final steps
    if (-e "$thread/result.fas") {                                                                            # Final steps
        system("cp", "$thread/result.fas", "$prefix\_result_$thread\.fas");                                   # Final steps
        # Print the name of the output file                                                                   # Final steps
        $print = "\tMatched region is saved to $folder/$prefix\_result_$thread\.fas\n";                       # Final steps
        print {$log} $print; print $print;                                                                    # Final steps
    }                                                                                                         # Final steps
                                                                                                              # Final steps
    # If the reference was matched to multiple contigs then copy the exonerate output file to the main folder # Final steps
    if ($multi{$thread} eq "Multiple contigs") {                                                              # Final steps
        # Worn the user that the matched region was not recovered and should be done manually                 # Final steps
        $print = "\tTry to recover the homologous region from the assembly\n";                                # Final steps
        print {$log} $print; print $print;                                                                    # Final steps
                                                                                                              # Final steps
        # Copy the exonerate output to the main folder and print its location                                 # Final steps
        system("cp", "$thread/$ref\.exonerate", "$prefix\_exonerate_$thread\.txt");                           # Final steps
        $print = "\tYou can find the exonerate output in $folder/$prefix\_exonerate_$thread\.txt\n";          # Final steps
        print {$log} $print; print $print;                                                                    # Final steps
    }                                                                                                         # Final steps
}                                                                                                             # Final steps
                                                                                                              # Final steps
# Delete every file except for the result files if hyperclean mode is selected                                # Final steps
system("rm", "-r", grep{$_ !~ /$prefix.*?\.((txt)|(fas)|(log))/} glob("*")) if $clean > 1;                    # Final steps
                                                                                                              # Final steps
# Print separator to mark the end of the run                                                                  # Final steps
$print = "=" x 90 . "\n" . "X" x 90 . "\n";                                                                   # Final steps
print {$log} $print; print $print;                                                                            # Final steps
###############################################################################################################


#=================================================SUBROUTINES================================================================================

###############################################################################################################
###############################################################################################################
##                                                SUBROUTINES                                                ##
###############################################################################################################
###############################################################################################################
# Defining variables                                                                                          #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # Defining variables
sub parse_input{                                                                                              # Defining variables
    # This subroutine reads the command line input: loads the variables, creates the folder and the files     # Defining variables
    #   Every input is a reference to a variable                                                              # Defining variables
    #   Outputs:    no return value                                                                           # Defining variables
    #               all the referenced variables maybe be modified, but the command line arguments            # Defining variables
    #   Created files:                                                                                        # Defining variables
    #                   1) the folder if needed                                                               # Defining variables
    #                   2) reference file                                                                     # Defining variables
    #                   3) bait file                                                                          # Defining variables
    #                   4) extra bait file if needed                                                          # Defining variables
    #                   5) symbolic links to the read files                                                   # Defining variables
    #   Inputs:                                                                                               # Defining variables
    my ($in_ref,        #  1) (array)  of arguments (@ARGV)                                                   # Defining variables
        $reference_ref, #  2) (scalar) to hold the location of the reference file                             # Defining variables
        $reads_ref,     #  3) (array)  to hold the location of the read files                                 # Defining variables
        $single_ref,    #  4) (scalar) which is TRUE if the reads are not paired                              # Defining variables
        $folder_ref,    #  5) (scalar) for the relative path of the folder to be used as main directory       # Defining variables
        $prefix_ref,    #  6) (scalar) to hold the prefix for the output files                                # Defining variables
        $type_ref,      #  7) (scalar) to hold the type of the run (defined by the --type argument)           # Defining variables
        $parameters_ref,#  8) (array)  to hold all the arguments that need to be handed to the assembler      # Defining variables
    #                               [0] name of the assembler program                                         # Defining variables
    #                               [1] additional parameters for the graph/hash generation                   # Defining variables
    #                               [2] additional parameters for the assembly generation                     # Defining variables
        $pre_log_ref,   #  9) (scalar) to hold all the test for the log file, until it is created             # Defining variables 
        $bait_ref,      # 10) (scalar) to hold the location of the bait file                                  # Defining variables
        $extra_bait_ref,# 11) (scalar) TRUE if there is was a bait defined in multi mode                      # Defining variables
        $clean_ref,     # 12) (scalar) defined if clean mode is selected,                                     # Defining variables
    #                                   it is 2 if hyperclean mode should be used                             # Defining variables
        $gzip_ref,      # 13) (scalar) TRUE if the read files are compressed                                  # Defining variables
        $min_len_ref,   # 14) (scalar) the minimal size for a contig to be kept                               # Defining variables
        $external_ref,  # 15) (scalar) the perl script to be used for the assembly                            # Defining variables
        $format_ref,    # 16) (scalar) the format of the read files                                           # Defining variables
                        ) = @_;                                                                               # Defining variables
                                                                                                              # Defining variables
    # Create a variable to hold the reference to the value that is waiting for a value                        # Defining variables
    my $current;                                                                                              # Defining variables
    # True if there should be a value (argument without "--" prefix)                                          # Defining variables
    my $wait;                                                                                                 # Defining variables
    # True if there could be more than one value                                                              # Defining variables
    my $array;                                                                                                # Defining variables
    # Variables to hold the information for the parameter array elements before they are stored in the array  # Defining variables
    my ($asmblr, $arg1, $arg2);  # @$parameters_ref = ($asmblr, $arg1, $arg2)                                 # Defining variables
                                                                                                              # Defining variables
    ########################################################################################################### Defining variables
    # Save the command line arguments into variables                                                          # Defining variables
    ########################################################################################################### Defining variables
    # Loop through the command line arguments and save the values to the correct place                        # Defining variables
    for (@$in_ref) {                                                                                          # Defining variables: get params
        if (/^--/) {                                                                                          # Defining variables: get params
            # If one of the arguments is waiting for a value then there is a problem with the input           # Defining variables: get params
            die "Parameter value is missing\n" if $wait;                                                      # Defining variables: get params
                                                                                                              # Defining variables: get params
            if (/--ref/) {                                                                                    # Defining variables: get params
                $current = $reference_ref;                                                                    # Defining variables: get params
                $array = undef;                                                                               # Defining variables: get params
                $wait++;                                                                                      # Defining variables: get params
                                                                                                              # Defining variables: get params
            } elsif (/--reads/) {                                                                             # Defining variables: get params
                $current = $reads_ref;                                                                        # Defining variables: get params
                $array++;                                                                                     # Defining variables: get params
                $wait++;                                                                                      # Defining variables: get params
                                                                                                              # Defining variables: get params
            } elsif (/--folder/) {                                                                            # Defining variables: get params
                $current = $folder_ref;                                                                       # Defining variables: get params
                $array = undef;                                                                               # Defining variables: get params
                $wait++;                                                                                      # Defining variables: get params
                                                                                                              # Defining variables: get params
            } elsif (/--prefix/) {                                                                            # Defining variables: get params
                $current = $prefix_ref;                                                                       # Defining variables: get params
                $array = undef;                                                                               # Defining variables: get params
                $wait++;                                                                                      # Defining variables: get params
                                                                                                              # Defining variables: get params
            } elsif (/--type/) {                                                                              # Defining variables: get params
                $current = $type_ref;                                                                         # Defining variables: get params
                $array = undef;                                                                               # Defining variables: get params
                $wait++;                                                                                      # Defining variables: get params
                                                                                                              # Defining variables: get params
            } elsif (/--arg1/) {                                                                              # Defining variables: get params
                $current = \$arg1;                                                                            # Defining variables: get params
                $array = undef;                                                                               # Defining variables: get params
                $wait++;                                                                                      # Defining variables: get params
                                                                                                              # Defining variables: get params
            } elsif (/--arg2/) {                                                                              # Defining variables: get params
                $current = \$arg2;                                                                            # Defining variables: get params
                $array = undef;                                                                               # Defining variables: get params
                $wait++;                                                                                      # Defining variables: get params
                                                                                                              # Defining variables: get params
            } elsif (/--assembler/) {                                                                         # Defining variables: get params
                $current = \$asmblr;                                                                          # Defining variables: get params
                $array = undef;                                                                               # Defining variables: get params
                $wait++;                                                                                      # Defining variables: get params
                                                                                                              # Defining variables: get params
            } elsif (/--single/) {                                                                            # Defining variables: get params
                $$single_ref++;                                                                               # Defining variables: get params
                $current = undef;                                                                             # Defining variables: get params
                $array = undef;                                                                               # Defining variables: get params
                                                                                                              # Defining variables: get params
            } elsif (/--clean/) {                                                                             # Defining variables: get params
                $$clean_ref++;                                                                                # Defining variables: get params
                $current = undef;                                                                             # Defining variables: get params
                $array = undef;                                                                               # Defining variables: get params
                                                                                                              # Defining variables: get params
            } elsif (/--bait/) {                                                                              # Defining variables: get params
                $current = $bait_ref;                                                                         # Defining variables: get params
                $array = undef;                                                                               # Defining variables: get params
                $$extra_bait_ref++;                                                                           # Defining variables: get params
                $wait++;                                                                                      # Defining variables: get params
                                                                                                              # Defining variables: get params
            } elsif (/--min_length=(\d+)/) {                                                                  # Defining variables: get params
                $$min_len_ref = $1;                                                                           # Defining variables: get params
                $current = undef;                                                                             # Defining variables: get params
                $array = undef;                                                                               # Defining variables: get params
                                                                                                              # Defining variables: get params
            } else {                                                                                          # Defining variables: get params
                # The given argument is not recognized so die                                                 # Defining variables: get params
                die "Unknown parameter encountered \"$_\"\n";                                                 # Defining variables: get params
            }                                                                                                 # Defining variables: get params
        } else {                                                                                              # Defining variables: get params
            # There was no argument initialized that is waiting for a value                                   # Defining variables: get params
            unless ($current) { die "Missing parameter value\n"}                                              # Defining variables: get params
                                                                                                              # Defining variables: get params
            # Store away the value according to its type (array or scalar)                                    # Defining variables: get params
            if ($array) {                                                                                     # Defining variables: get params
                push @$current, $_;                                                                           # Defining variables: get params
            } else {                                                                                          # Defining variables: get params
                $$current = $_;                                                                               # Defining variables: get params
            }                                                                                                 # Defining variables: get params
                                                                                                              # Defining variables: get params
            # Except for array variables there is no more value needed                                        # Defining variables: get params
            $wait = undef;                                                                                    # Defining variables: get params
            # if the assembler is external than wait for the next word                                        # Defining variables: get params
            if ($_ eq "external") {                                                                           # Defining variables: get params
                $current = $external_ref;                                                                     # Defining variables: get params
                $wait++;                                                                                      # Defining variables: get params
            }                                                                                                 # Defining variables: get params
        }                                                                                                     # Defining variables: get params
    ########################################################################################################### Defining variables: get params
    # End of the for loop of the command line arguments                                                       # Defining variables: get params
    ########################################################################################################### Defining variables: get params
    }                                                                                                         # Defining variables: get params
                                                                                                              # Defining variables
    # Check if all the necessary variables are defined and the files are reachable                            # Defining variables: essentials
    #   reference file                                                                                        # Defining variables: essentials
    unless ($$reference_ref) {                                                                                # Defining variables: essentials
        print "\n$usage\n";                                                                                   # Defining variables: essentials
        die "Necessary argument is missing: no reference file was defined\n";                                 # Defining variables: essentials
    }                                                                                                         # Defining variables: essentials
    # Check whether the specified file exists                                                                 # Defining variables: essentials
    if (-e $$reference_ref) {                                                                                 # Defining variables: essentials
        # Print the name of the reference file to the log and STDOUT                                          # Defining variables: essentials
        $$pre_log_ref .= "The reference file is $$reference_ref\n";                                           # Defining variables: essentials
    } else {                                                                                                  # Defining variables: essentials
        die "The specified reference file ($$reference_ref) is missing\n";                                    # Defining variables: essentials
    }                                                                                                         # Defining variables: essentials
                                                                                                              # Defining variables: essentials
    #   read files                                                                                            # Defining variables: essentials
    if (scalar @$reads_ref == 0) {                                                                            # Defining variables: essentials
        print "\n$usage\n";                                                                                   # Defining variables: essentials
        die "Necessary argument is missing: no read file was defined\n";                                      # Defining variables: essentials
    }                                                                                                         # Defining variables: essentials
    # Check whether the specified files exist                                                                 # Defining variables: essentials
    for (@$reads_ref) {                                                                                       # Defining variables: essentials
        if (-e $_) {                                                                                          # Defining variables: essentials
            # Print the name of the read file to the log and STDOUT                                           # Defining variables: essentials
            $$pre_log_ref .= "The read file is $_\n";                                                         # Defining variables: essentials
        } else {                                                                                              # Defining variables: essentials
            die "The specified read file ($_) is missing\n";                                                  # Defining variables: essentials
        }                                                                                                     # Defining variables: essentials
    }                                                                                                         # Defining variables: essentials
                                                                                                              # Defining variables: essentials
    #   folder                                                                                                # Defining variables: essentials
    unless ($$folder_ref) {                                                                                   # Defining variables: essentials
        print "\n$usage\n";                                                                                   # Defining variables: essentials
        die "Necessary argument is missing: no folder was defined\n";                                         # Defining variables: essentials
    }                                                                                                         # Defining variables: essentials
                                                                                                              # Defining variables: essentials
    #   prefix                                                                                                # Defining variables: essentials
    unless ($$prefix_ref) {                                                                                   # Defining variables: essentials
        print "\n$usage\n";                                                                                   # Defining variables: essentials
        die "Necessary argument is missing: no prefix was defined\n";                                         # Defining variables: essentials
    }                                                                                                         # Defining variables: essentials
                                                                                                              # Defining variables
    # Extra bait is only important to know if in multi mode, because then there should be an inner loop       # Defining variables: extra bait
    if ($$type_ref) {                                                                                         # Defining variables: extra bait
        if ($$type_ref !~ /multi/) {                                                                          # Defining variables: extra bait
            $$extra_bait_ref = undef;                                                                         # Defining variables: extra bait
        }                                                                                                     # Defining variables: extra bait
    } else {                                                                                                  # Defining variables: extra bait
        $$extra_bait_ref = undef;                                                                             # Defining variables: extra bait
    }                                                                                                         # Defining variables: extra bait
                                                                                                              # Defining variables
    # If an extra bait file was defined then combine it with the reference to use as bait                     # Defining variables: Bait
    #   ELSE use reference as bait                                                                            # Defining variables: Bait
    if ($$bait_ref) {                                                                                         # Defining variables: Bait
        # Check if the file really exists or not                                                              # Defining variables: Bait
        if (-e $$bait_ref) {                                                                                  # Defining variables: Bait
            # Print the name of the file to the log and STDOUT                                                # Defining variables: Bait
            $$pre_log_ref .= "The first bait is a concatenation of $$bait_ref and $$reference_ref\n";         # Defining variables: Bait
        } else {                                                                                              # Defining variables: Bait
            die "The specified bait file ($$bait_ref) is missing\n";                                          # Defining variables: Bait
        }                                                                                                     # Defining variables: Bait
    } else {                                                                                                  # Defining variables: Bait
        # Use the reference file as bait                                                                      # Defining variables: Bait
        $$bait_ref = $$reference_ref;                                                                         # Defining variables: Bait
        $$pre_log_ref .= "The first bait file is $$reference_ref\n";                                          # Defining variables: Bait
    }                                                                                                         # Defining variables
                                                                                                              # Defining variables
    # Check if single or paired mode should be used                                                           # Defining variables: Single mode
    #   Paired mode if there are two read files and --single was not declared                                 # Defining variables: Single mode
    if (scalar @$reads_ref == 2 && not $$single_ref) {                                                        # Defining variables: Single mode
        $$pre_log_ref .= "Using paired-end mode\n";                                                           # Defining variables: Single mode
    } else {                                                                                                  # Defining variables: Single mode
        $$single_ref++;                                                                                       # Defining variables: Single mode
        $$pre_log_ref .= "Using single-end mode\n";                                                           # Defining variables: Single mode
    }                                                                                                         # Defining variables
                                                                                                              # Defining variables
    # Sort out the assembler parameters                                                                       # Defining variables: assembler
    #   If no assembler was specified then use edena                                                          # Defining variables: assembler
    $asmblr = "edena" unless $asmblr;                                                                         # Defining variables: assembler
    #   Check if the assembler specified is an acceptable value or not                                        # Defining variables: assembler
    unless ($asmblr eq "edena" || $asmblr eq "velvet" || $asmblr eq "alternative" || $asmblr eq "external") { # Defining variables: assembler
        die "Unknown assembler specified, possible values are: " .                                            # Defining variables: assembler
            "'edena', 'velvet', 'alternative' or 'external'\n";                                               # Defining variables: assembler
    }                                                                                                         # Defining variables: assembler
                                                                                                              # Defining variables: assembler
    # Put the assembler values in the parameter array                                                         # Defining variables: assembler
    @$parameters_ref = ($asmblr, $arg1, $arg2);                                                               # Defining variables: assembler
                                                                                                              # Defining variables
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - # Defining variables: Folder
    # Check the folder                                                                                        # Defining variables: Folder
    #   IF it is not present it needs to be created                                                           # Defining variables: Folder
    #   ELSE check if empty                                                                                   # Defining variables: Folder
    #       ELSE ask if you want    1) to use the files in there                                              # Defining variables: Folder
    #                               2) to delete the content                                                  # Defining variables: Folder
    #                               3) to quit                                                                # Defining variables: Folder
    # Remove trailing "/" form the folder name                                                                # Defining variables: Folder
    $$folder_ref =~ s/\/$//;                                                                                  # Defining variables: Folder
    if (-d $$folder_ref) {                                                                                    # Defining variables: Folder
        # Get the content of the folder                                                                       # Defining variables: Folder
        my @list = glob("$$folder_ref/*");                                                                    # Defining variables: Folder
        # If there is something in it then ask what to do                                                     # Defining variables: Folder
        if (@list) {                                                                                          # Defining variables: Folder
            print "Specified folder ($$folder_ref) is not empty.\n";                                          # Defining variables: Folder
            print "Do you wish to continue?[y/n] ";                                                           # Defining variables: Folder
            # Get the response of the user                                                                    # Defining variables: Folder
            my $response = <STDIN>;                                                                           # Defining variables: Folder
            # If it is yes then continue else quit the program                                                # Defining variables: Folder
            if ($response =~ /^y/i) {                                                                         # Defining variables: Folder
                # Ask whether the content should be reused or should it be deleted                            # Defining variables: Folder
                # If the content is not overlapping with the file naming of the program then it won't         # Defining variables: Folder
                #   interfere with the run. Beware: running hyperclean mode will delete the old content of    # Defining variables: Folder
                #   the folder at the end of the run                                                          # Defining variables: Folder
                print "Do you want the program to use the content of the folder?\n(Already existing files",   # Defining variables: Folder
                " will not be overwritten, this includes reference, bait and read files. ",                   # Defining variables: Folder
                "If you wish to replace them delete those files, then run the command again) [y/n] ";         # Defining variables: Folder
                # Get the user's response                                                                     # Defining variables: Folder
                $response = <STDIN>;                                                                          # Defining variables: Folder
                if ($response =~ /^y/i) {                                                                     # Defining variables: Folder
                    # Print to the log that the selected folder is chosen and that it is not-empty            # Defining variables: Folder
                    $$pre_log_ref .= "The working directory is $$folder_ref (the folder is not empty)\n";     # Defining variables: Folder
                } else {                                                                                      # Defining variables: Folder
                    # Ask if the content should be deleted or not                                             # Defining variables: Folder
                    print "Do you wish to delete the content of the folder?[y/n] ";                           # Defining variables: Folder
                    $response = <STDIN>;                                                                      # Defining variables: Folder
                    if ($response =~ /^y/i) {                                                                 # Defining variables: Folder
                        # Just to make sure ask it for a second time                                          # Defining variables: Folder
                        print "Are you sure?[y/n] ";                                                          # Defining variables: Folder
                        $response = <STDIN>;                                                                  # Defining variables: Folder
                        if ($response =~ /^y/i) {                                                             # Defining variables: Folder
                            # Delete the folder and its content                                               # Defining variables: Folder
                            system("rm", "-r", $$folder_ref);                                                 # Defining variables: Folder
                            # Create the folder anew                                                          # Defining variables: Folder
                            mkdir $$folder_ref;                                                               # Defining variables: Folder
                            # Print to the log that this folder is used and it was emptied                    # Defining variables: Folder
                            $$pre_log_ref .= "The working directory is $$folder_ref (content deleted)\n";     # Defining variables: Folder
                        } else {                                                                              # Defining variables: Folder
                            die "The program is terminated\n";                                                # Defining variables: Folder
                        }                                                                                     # Defining variables: Folder
                    } else {                                                                                  # Defining variables: Folder
                        die "The program is terminated\n";                                                    # Defining variables: Folder
                    }                                                                                         # Defining variables: Folder
                }                                                                                             # Defining variables: Folder
            } else {                                                                                          # Defining variables: Folder
                die "The program is terminated\n";                                                            # Defining variables: Folder
            }                                                                                                 # Defining variables: Folder
        }                                                                                                     # Defining variables: Folder
    } else {                                                                                                  # Defining variables: Folder
        # Create folder and print the name to the log                                                         # Defining variables: Folder
        mkdir $$folder_ref;                                                                                   # Defining variables: Folder
        $$pre_log_ref .= "The working directory is $$folder_ref\n";                                           # Defining variables: Folder
    }                                                                                                         # Defining variables: Folder
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - # Defining variables: Folder
    # Create a copy of the perl script of the assembler if it was specified                                   # Defining variables: Copy assembler
    if ($asmblr eq "external") {                                                                              # Defining variables: Copy assembler
        system("cp", $$external_ref, "$$folder_ref\/assembler.pl");                                           # Defining variables: Copy assembler
    }                                                                                                         # Defining variables: Copy assembler
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - # Defining variables: Create
    # Create copies and symlinks inside the working directory to be used by the program                       # Defining variables: Create
    #   Create the reference file if it isn't present already                                                 # Defining variables: Create ref
    unless (-e "$$folder_ref\/reference.fas") {                                                               # Defining variables: Create ref
        # Copy the specified reference into the main folder                                                   # Defining variables: Create ref
        system("cp", $$reference_ref, "$$folder_ref\/reference.fas");                                         # Defining variables: Create ref
        # Print that the reference was copied into the folder                                                 # Defining variables: Create ref
        $$pre_log_ref .= "\tReference file is copied to $$folder_ref\/reference.fas\n";                       # Defining variables: Create ref
    } else {                                                                                                  # Defining variables: Create ref
        # Print that the previous reference file (found in the folder) is used and not the one defined        # Defining variables: Create ref
        #   at invocation                                                                                     # Defining variables: Create ref
        $$pre_log_ref .= "\tPreviously defined reference file is used $$folder_ref\/reference.fas\n";         # Defining variables: Create ref
    }                                                                                                         # Defining variables: Create ref
                                                                                                              # Defining variables: Create
    # Create the bait file if it is not present yet                                                           # Defining variables: Create bait
    unless (-e "$$folder_ref\/bait.fas") {                                                                    # Defining variables: Create bait
        # Combine the reference the bait file (if specified) together to create the first bait file           # Defining variables: Create bait
        my @baits = ($$bait_ref, $$reference_ref);                                                            # Defining variables: Create bait
        `cat @baits | perl -ne 'if (/>/) {\$n++; \$_ = ">\$n\n"}; print' >$$folder_ref\/bait.fas `;           # Defining variables: Create bait
                                                                                                              # Defining variables: Create bait
        # copy the bait file into the folder so it is preserved within the folder (will not be used)          # Defining variables: Create bait
        system("cp", $$bait_ref, "$$folder_ref\/extra_bait.fas");                                             # Defining variables: Create bait
        # Print that the bait was create inside the folder                                                    # Defining variables: Create bait
        $$pre_log_ref .= "\tBait file is created: $$folder_ref\/bait.fas\n";                                  # Defining variables: Create bait
    } else {                                                                                                  # Defining variables: Create bait
        # Print that the previous reference file (found in the folder) is used and not the one defined        # Defining variables: Create bait
        #   at invocation                                                                                     # Defining variables: Create bait
        $$pre_log_ref .= "\tPreviously defined  bait file is used together with the new files\n";             # Defining variables: Create bait
        my @baits = ($$bait_ref, $$reference_ref, "$$folder_ref\/bait.fas");                                  # Defining variables: Create bait
        `cat @baits | perl -ne 'if (/>/) {\$n++; \$_ = ">\$n\n"}; print' >$$folder_ref\/temp_bait.fas `;      # Defining variables: Create bait
        system('cp', "$$folder_ref\/temp_bait.fas", "$$folder_ref\/bait.fas");                                # Defining variables: Create bait
    }                                                                                                         # Defining variables: Create bait
                                                                                                              # Defining variables
    # Set the names of the bait and reference file                                                            # Defining variables
    $$bait_ref = "bait.fas";                                                                                  # Defining variables
    $$reference_ref = "reference.fas";                                                                        # Defining variables
                                                                                                              # Defining variables
    # Create symlinks to the read files and store the names of the links in an array                          # Defining variables: Create reads
    #   This enables us to use the same array to refer to the reads in all the iterations                     # Defining variables: Create reads
    # Initialize a counter for the read files. This will keep track of what is the serial number of the       # Defining variables: Create reads
    #   current read file                                                                                     # Defining variables: Create reads
    my $i;                                                                                                    # Defining variables: Create reads
    for (@$reads_ref) {                                                                                       # Defining variables: Create reads
        $i++;                                                                                                 # Defining variables: Create reads
        # If the specified read file is zipped then increment the gzip_ref value                              # Defining variables: Create reads
        $$gzip_ref++ if /\.gz/;                                                                               # Defining variables: Create reads
                                                                                                              # Defining variables: Create reads
        # Identify the file format of the read file based on the first character                              # Defining variables: Create reads
        my $read_cmd = "cat";                                                                                 # Defining variables: Create reads
        $read_cmd = "zcat" if $$gzip_ref;                                                                     # Defining variables: Create reads
        my $first_line = `$read_cmd $_ | head -1`;                                                            # Defining variables: Create reads
        if ($first_line =~ /^>/) {                                                                            # Defining variables: Create reads
            if ($$format_ref) {                                                                               # Defining variables: Create reads
                die "All the read files should have the same format\n" if $$format_ref ne "fasta";            # Defining variables: Create reads
            }                                                                                                 # Defining variables: Create reads
            $$format_ref = "fasta";                                                                           # Defining variables: Create reads
        } elsif ($first_line =~ /^\@/) {                                                                      # Defining variables: Create reads
            if ($$format_ref) {                                                                               # Defining variables: Create reads
                die "All the read files should have the same format\n" if $$format_ref ne "fastq";            # Defining variables: Create reads
            }                                                                                                 # Defining variables: Create reads
            $$format_ref = "fastq";                                                                           # Defining variables: Create reads
        } else {                                                                                              # Defining variables: Create reads
                                                                                                              # Defining variables: Create reads
        }                                                                                                     # Defining variables: Create reads
                                                                                                              # Defining variables: Create reads
        # Create the name of the symlink for the read files                                                   # Defining variables: Create reads
        my $symlink = "$$folder_ref\/reads$i\.$$format_ref";                                                  # Defining variables: Create reads
                                                                                                              # Defining variables: Create reads
        if ($$gzip_ref) {                                                                                     # Defining variables: Create reads
            # If there is a zipped read file then all the read files should be zipped files                   # Defining variables: Create reads
            die "All the reads should be either gzipped or uncompressed\n" unless $$gzip_ref == $i;           # Defining variables: Create reads
        }                                                                                                     # Defining variables: Create reads
                                                                                                              # Defining variables: Create reads
        # Create a symlink unless there is already on inside the folder                                       # Defining variables: Create reads
        unless (-l $symlink) {                                                                                # Defining variables: Create reads
            # Create the symlink to the specified read files                                                  # Defining variables: Create reads
            symlink(File::Spec->rel2abs($_), $symlink);                                                       # Defining variables: Create reads
                                                                                                              # Defining variables: Create reads
            # Check if the link is present and that it is working or die                                      # Defining variables: Create reads
            if (-l $symlink && stat($symlink)) {                                                              # Defining variables: Create reads
                # Print the name of the symlink created                                                       # Defining variables: Create reads
                $$pre_log_ref .= "\tSymlink ($symlink) is created for the read file ($_)\n";                  # Defining variables: Create reads
            } else {                                                                                          # Defining variables: Create reads
                die "Could not create symlinks for a read file ($_).\nProgram is terminated\n";               # Defining variables: Create reads
            }                                                                                                 # Defining variables: Create reads
        } else {                                                                                              # Defining variables: Create reads
            # Print that the previously created symlink will be used                                          # Defining variables: Create reads
            $$pre_log_ref .= "\tPreviously created symlink is used ($$folder_ref\/reads$i\.$$format_ref)\n";  # Defining variables: Create reads
        }                                                                                                     # Defining variables: Create reads
                                                                                                              # Defining variables: Create reads
        # Change the relative path of the symlink that will be used during the program to the read file array # Defining variables: Create reads
        $$reads_ref[$i - 1] = "../reads$i\.$$format_ref";                                                     # Defining variables: Create reads
    }                                                                                                         # Defining variables: Create reads
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - # Defining variables
}                                                                                                             # Defining variables
###############################################################################################################
#############################################################################
# Update the path of variables after changing directories                   #
#   First input is the prefix for the files                                 #
#   This subroutine takes the reference of the following files:             #
#       reference, bait, old collection                                     #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
sub update_path{                                                            #
    my ($path, $ref_ref, $bait_ref, $old_ref) = @_;                         #
    $$ref_ref =  $$path . "reference.fas";                                  #
    $$bait_ref = $$path . "bait.fas";                                       #
    $$old_ref = $$path . "old_collection.list";                             #
}                                                                           #
#############################################################################
#####################################################################################################
# Subroutines to handle the multi hash                                                              # Multi thread
#   %multi:  keys:   the id number of the given thread/run                                          # Multi thread
#            values: if active => FALSE else TRUE (a string describing the why it has terminated)   # Multi thread
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
sub load_multi{                                                                                     # Multi thread (load_multi)
    # This subroutine creates the keys and values for the multi hash variable                       # Multi thread (load_multi)
    #   If multi mode was specified, then the reference will be split up to each entry and a single # Multi thread (load_multi)
    #       thread will be created for it in the multi hash                                         # Multi thread (load_multi)
    #   Also it creates the folders for each of the threads in the main folder and creats the       # Multi thread (load_multi)
    #       reference, bait files (they are the same) and old collection for each thread, unless    # Multi thread (load_multi)
    #       they already exist                                                                      # Multi thread (load_multi)
    #   Inputs:                                                                                     # Multi thread (load_multi)
    my ($log_ref,   # The reference of the log handle                                               # Multi thread (load_multi)
        $ref_ref,   # The reference of the location of the reference file                           # Multi thread (load_multi)
        $bait_ref,  # The reference of the location of the bait file                                # Multi thread (load_multi)
	$type_ref,  # The reference of the run type variable                                        # Multi thread (load_multi)
        $multi_ref, # The reference to the multi hash variable                                      # Multi thread (load_multi)
        $types_ref, # The reference to the types hash variable                                      # Multi thread (load_multi)
                    ) = @_;                                                                         # Multi thread (load_multi)
                                                                                                    # Multi thread (load_multi)
    # A hash to store the entries of the reference file as key-value pairs                          # Multi thread (load_multi)
    #   key: id and value: sequence                                                                 # Multi thread (load_multi)
    my %ref_fasta;                                                                                  # Multi thread (load_multi)
    # An array to store the proper order of the entries in the reference file                       # Multi thread (load_multi)
    my @order;                                                                                      # Multi thread (load_multi)
                                                                                                    # Multi thread (load_multi)
    # If multi type is specified then create a separate thread for each of the entries in the       # Multi thread (load_multi)
    #   reference file                                                                              # Multi thread (load_multi)
    if ($$type_ref && $$type_ref =~ /multi/) {                                                      # Multi thread (load_multi)
        # Print that multi thread mode is selected                                                  # Multi thread (load_multi)
        $print = "Multi thread mode is selected: Reference file is going to be split up\n";         # Multi thread (load_multi)
        print {$$log_ref} $print; print $print;                                                     # Multi thread (load_multi)
                                                                                                    # Multi thread (load_multi)
        # Create folders in the main folder for each reference entry                                # Multi thread (load_multi)
        # Create the necessary files within that folder                                             # Multi thread (load_multi)
        #   old_collection.list     -empty file, will be updated later                              # Multi thread (load_multi)
        #   reference.fas           -one of the fasta entries in the reference file                 # Multi thread (load_multi)
        #   bait.fas                -same as the reference.fas, but will be updated later           # Multi thread (load_multi)
                                                                                                    # Multi thread (load_multi)
        # Check if the reference file is correct (no entry with the same id)                        # Multi thread (load_multi)
        &check_fasta($log_ref, $ref_ref);                                                           # Multi thread (load_multi)
                                                                                                    # Multi thread (load_multi)
        # Load the content of the reference file into a hash and store the order of the entries     # Multi thread (load_multi)
        &fasta2hash($ref_ref, \%ref_fasta, \@order);                                                # Multi thread (load_multi)
                                                                                                    # Multi thread (load_multi)
        my $i;                                                                                      # Multi thread (load_multi)
        for (@order) {                                                                              # Multi thread (load_multi)
            $i++;                                                                                   # Multi thread (load_multi)
            my $dir = "thread_$i";                                                                  # Multi thread (load_multi)
            if ($_ =~ /(exonerate)/) {                                                              # Multi thread (load_multi)
                $$types_ref{$dir} = $1;                                                             # Multi thread (load_multi)
            }                                                                                       # Multi thread (load_multi)
            if ($_ =~ /(total=\d+)/) {                                                              # Multi thread (load_multi)
                $$types_ref{$dir} .= $1;                                                            # Multi thread (load_multi)
            }                                                                                       # Multi thread (load_multi)
            if ($_ =~ /(longest=\d+)/) {                                                            # Multi thread (load_multi)
                $$types_ref{$dir} .= $1;                                                            # Multi thread (load_multi)
            }                                                                                       # Multi thread (load_multi)
            if ($_ =~ /(N50=\d+)/) {                                                                # Multi thread (load_multi)
                $$types_ref{$dir} .= $1;                                                            # Multi thread (load_multi)
            }                                                                                       # Multi thread (load_multi)
            if ($_ =~ /(exhaustive)/) {                                                             # Multi thread (load_multi)
                $$types_ref{$dir} = $1;                                                             # Multi thread (load_multi)
            }                                                                                       # Multi thread (load_multi)
            $$types_ref{$dir} = $$type_ref unless $$types_ref{$dir};                                # Multi thread (load_multi)
            # create the hash entry, value is undefined => active thread                            # Multi thread (load_multi)
            $$multi_ref{$dir} = undef;                                                              # Multi thread (load_multi)
            mkdir $dir unless -d $dir;                                                              # Multi thread (load_multi)
            open my $empty, '>>', "$dir\/old_collection.list";                                      # Multi thread (load_multi)
            close $empty;                                                                           # Multi thread (load_multi)
            if (-e "$dir\/reference.fas") {                                                         # Multi thread (load_multi)
                $print = "\tThe reference file of $dir is replaced\n";                              # Multi thread (load_multi)
                print {$$log_ref} $print; print $print;                                             # Multi thread (load_multi)
            } else {                                                                                # Multi thread (load_multi)
                $print = "\tThe reference file of $dir is created\n";                               # Multi thread (load_multi)
                print {$$log_ref} $print; print $print;                                             # Multi thread (load_multi)
            }                                                                                       # Multi thread (load_multi)
            $print = "\t\t$_\n";                                                                    # Multi thread (load_multi)
            print {$$log_ref} $print; print $print;                                                 # Multi thread (load_multi)
            open my $ref_handle, '>', "$dir\/reference.fas";                                        # Multi thread (load_multi)
            print {$ref_handle} &to_fasta($_, $ref_fasta{$_});                                      # Multi thread (load_multi)
            close $ref_handle;                                                                      # Multi thread (load_multi)
            `cat $dir\/reference.fas >>$dir\/bait.fas`;                                             # Multi thread (load_multi)
        }                                                                                           # Multi thread (load_multi)
    } else {                                                                                        # Multi thread (load_multi)
        $print = "Single thread mode is selected: Reference file is kept intact\n";                 # Multi thread (load_multi)
        print {$$log_ref} $print; print $print;                                                     # Multi thread (load_multi)
        my $dir = "thread_1";                                                                       # Multi thread (load_multi)
        # create the hash entry, value is undefined => active thread                                # Multi thread (load_multi)
        $$multi_ref{$dir} = undef;                                                                  # Multi thread (load_multi)
        mkdir $dir unless -d $dir;                                                                  # Multi thread (load_multi)
        open my $empty, '>>', "$dir\/old_collection.list";                                          # Multi thread (load_multi)
        close $empty;                                                                               # Multi thread (load_multi)
        if (-e "$dir\/reference.fas") {                                                             # Multi thread (load_multi)
            $print = "\tThe reference file of $dir is replaced\n";                                  # Multi thread (load_multi)
            print {$$log_ref} $print; print $print;                                                 # Multi thread (load_multi)
        } else {                                                                                    # Multi thread (load_multi)
            $print = "\tThe reference file of $dir is created\n";                                   # Multi thread (load_multi)
            print {$$log_ref} $print; print $print;                                                 # Multi thread (load_multi)
        }                                                                                           # Multi thread (load_multi)
        system("cp", $$ref_ref, "$dir\/reference.fas");                                             # Multi thread (load_multi)
        `cat $$bait_ref >>$dir\/bait.fas`;                                                          # Multi thread (load_multi)
    }                                                                                               # Multi thread (load_multi)
}                                                                                                   # Multi thread (load_multi)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
sub get_actives{                                                                                    # Multi thread (get_actives)
    # This subroutine has only one input (the multi hash)                                           # Multi thread (get_actives)
    #	and returns the active threads (with undefined values) in a sorted array                    # Multi thread (get_actives)
    my ($hash_ref) = @_;                                                                            # Multi thread (get_actives)
                                                                                                    # Multi thread (get_actives)
    # Get all the name of the active threads from the hash                                          # Multi thread (get_actives)
    my @unsorted = grep{not $$hash_ref{$_}} keys %$hash_ref;                                        # Multi thread (get_actives)
                                                                                                    # Multi thread (get_actives)
    # If there is more then one active thread then sort them by their number                        # Multi thread (get_actives)
    if (1 < scalar @unsorted) {                                                                     # Multi thread (get_actives)
        my @sorted = sort{  my ($a1) = $a =~ /thread_(\d+)/;    # get the number of the thread      # Multi thread (get_actives)
                            my ($b1) = $b =~ /thread_(\d+)/;    # get the number of the thread      # Multi thread (get_actives)
                            $a1 <=> $b1                         } @unsorted;                        # Multi thread (get_actives)
        # Return the list                                                                           # Multi thread (get_actives)
        return @sorted;                                                                             # Multi thread (get_actives)
    } else {                                                                                        # Multi thread (get_actives)
        return @unsorted;                                                                           # Multi thread (get_actives)
    }                                                                                               # Multi thread (get_actives)
}                                                                                                   # Multi thread (get_actives)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
sub to_fasta {                                                                                      # Multi thread (to_fasta)
    # This subroutine takes two inputs: sequence id and sequence                                    # Multi thread (to_fasta)
    #    and returns a fasta formatted string                                                       # Multi thread (to_fasta)
    my ($seq_name, $seq) = @_;                                                                      # Multi thread (to_fasta)
    # The first line is the id line                                                                 # Multi thread (to_fasta)
    my $formatted_seq = ">$seq_name\n";                                                             # Multi thread (to_fasta)
    # The rest is the sequence, with a maximal line length of 60                                    # Multi thread (to_fasta)
    while (my $chunk = substr($seq, 0, 60, "")) {                                                   # Multi thread (to_fasta)
        $formatted_seq .= "$chunk\n";                                                               # Multi thread (to_fasta)
    }                                                                                               # Multi thread (to_fasta)
    # Return the string                                                                             # Multi thread (to_fasta)
    return $formatted_seq;                                                                          # Multi thread (to_fasta)
}                                                                                                   # Multi thread (to_fasta)
#####################################################################################################
########################################################################################################################
# Check all the necessary files if they are present, non-empty and correct                                             #
## # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # Check files
sub check_files{                                                                                                       # Check files
    # This subroutine checks if all the files are available and in proper format                                       # Check files
    #   Inputs:                                                                                                        # Check files
    my ($log_ref,       # The reference of the log file handle                                                         # Check files
        $bait_ref,      # The reference of the bait file                                                               # Check files
        $ref_ref,       # The reference of the reference file                                                          # Check files
        $reads_ref,     # The reference of the array storing the read files                                            # Check files
        $single_ref,    # The reference of the variable that is true if the reads are unpaired                         # Check files
        $gzip_ref       # The reference of the variable that is true if the read files are compressed by gzip          # Check files
                        ) = @_;                                                                                        # Check files
                                                                                                                       # Check files
    # Check the two fasta format files                                                                                 # Check files
    for ($$bait_ref, $$ref_ref) {                                                                                      # Check files
        # Check if they exist and if they are non-empty                                                                # Check files
        if (-e $_ && not -z $_) {                                                                                      # Check files
            # This function tests fasta file, that all the entries have unique identifiers                             # Check files
            &check_fasta($log_ref, \$_);                                                                               # Check files
        } else {                                                                                                       # Check files
            # Print an error if they are missing or are empty                                                          # Check files
            $print = "Error: the fasta file ($_) does not exist or is empty!\n";                                       # Check files
            print {$$log_ref} $print; print $print;                                                                    # Check files
            die "Fasta file ($_) does not exist or is empty\nProgram is terminated\n";                                 # Check files
        }                                                                                                              # Check files
    }                                                                                                                  # Check files
                                                                                                                       # Check files
    # Check the read file symlinks if they are working properly                                                        # Check files
    for my $read (@$reads_ref) {                                                                                       # Check files
        # The symlink exists and is working                                                                            # Check files
        if (-l $read && stat($read)) {                                                                                 # Check files
            # The file it is linked to is non-empty                                                                    # Check files
            if (-z readlink($read)) {                                                                                  # Check files
                # Print an error if it is empty                                                                        # Check files
                $print = "Error: the symlink ($read) points to an empty file(" . readlink($read) . ")\n";              # Check files
                print {$$log_ref} $print; print $print;                                                                # Check files
                die "Symlink ($read) points to an empty file (" . readlink($read) . ")\nProgram is terminated\n";      # Check files
            }                                                                                                          # Check files
        } else {                                                                                                       # Check files
            # Print an error if the symlink is missing or it isn't working                                             # Check files
            $print = "Error: the symlink ($read) to one of the read file does not exist or it is broken\n";            # Check files
            print {$$log_ref} $print; print $print;                                                                    # Check files
            die "Symlink ($read) does not exist or it is broken\nProgram is terminated\n";                             # Check files
        }                                                                                                              # Check files
    }                                                                                                                  # Check files
                                                                                                                       # Check files
    # Check the read files                                                                                             # Check files
    #	That they are readable, are in either fastq or fasta format and if paired then the pair can be recovered       # Check files
    &check_read_file($log_ref, $reads_ref, $single_ref, $gzip_ref);                                                    # Check files
}                                                                                                                      # Check files
## # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
sub check_fasta{                                                                                                       # Check files (fasta)
    # Check a fasta file for errors: duplicate entries, empty entries or incorrect characters                          # Check files (fasta)
    #   Inputs: reference for the log file handle and the reference of the fasta file                                  # Check files (fasta)
    my ($log_ref, $file_ref) = @_;                                                                                     # Check files (fasta)
                                                                                                                       # Check files (fasta)
    # Hash to store fasta data: key: first word of the id line, value: the sequence                                    # Check files (fasta)
    my %hash;                                                                                                          # Check files (fasta)
    # A scalar to store the sequence id when loading the sequence                                                      # Check files (fasta)
    my $cur;                                                                                                           # Check files (fasta)
                                                                                                                       # Check files (fasta)
    # Open a file handle for the fasta file                                                                            # Check files (fasta)
    open (my $fasta, '<', $$file_ref);                                                                                 # Check files (fasta)
    # Read in the fasta data and check for entries with same id line (only looking at the first words)                 # Check files (fasta)
    for(<$fasta>) {                                                                                                    # Check files (fasta)
        # If id line get the first word, because mirabait and exonerate use the first word only                        # Check files (fasta)
        if (/^>(\S+)/) {                                                                                               # Check files (fasta)
            # Save the id of the current entry                                                                         # Check files (fasta)
            $cur = $1;                                                                                                 # Check files (fasta)
                                                                                                                       # Check files (fasta)
            # Check if there is already a an entry in the hash with the same id                                        # Check files (fasta)
            if ($hash{$cur}) {                                                                                         # Check files (fasta)
                # Print an error and die if there is duplicate entry                                                   # Check files (fasta)
                $print = "Error: the fasta file ($$file_ref) has duplicate entry names ($cur)!\n";                     # Check files (fasta)
                print {$$log_ref} $print; print $print;                                                                # Check files (fasta)
                die "Fasta file ($$file_ref) has duplicate entry names!\nProgram is terminated\n";                     # Check files (fasta)
            }                                                                                                          # Check files (fasta)
        } else {                                                                                                       # Check files (fasta)
            # This a sequence line belonging to the id defined by the last id line                                     # Check files (fasta)
            # Remove white space and new line endings, so only the actual sequence is stored away                      # Check files (fasta)
            s/\s//g;                                                                                                   # Check files (fasta)
            s/\R//g;                                                                                                   # Check files (fasta)
                                                                                                                       # Check files (fasta)
            # Make sure that there are no incorrect characters in the sequences                                        # Check files (fasta)
            if (/[^ACTGRYMKSWHBVDN]/i) {                                                                               # Check files (fasta)
                $print = "Error: the fasta file ($$file_ref) has entry containing incorrect character!\n";             # Check files (fasta)
                print {$$log_ref} $print; print $print;                                                                # Check files (fasta)
                die "Fasta file ($$file_ref) contains incorrect character!\nProgram is terminated\n";                  # Check files (fasta)
            }                                                                                                          # Check files (fasta)
            # Append the sequence to the sequence of the current id                                                    # Check files (fasta)
            $hash{$cur} .= $_;                                                                                         # Check files (fasta)
        }                                                                                                              # Check files (fasta)
    }                                                                                                                  # Check files (fasta)
    # Close the fasta file handle                                                                                      # Check files (fasta)
    close $fasta;                                                                                                      # Check files (fasta)
                                                                                                                       # Check files (fasta)
    # Check for entries without a sequence                                                                             # Check files (fasta)
    for (keys %hash) {                                                                                                 # Check files (fasta)
        unless ($hash{$_}) {                                                                                           # Check files (fasta)
            # Print an error and die, because there is an id without sequence                                          # Check files (fasta)
            $print = "Error: the fasta file ($$file_ref) has empty entry ($_)!\n";                                     # Check files (fasta)
            print {$$log_ref} $print; print $print;                                                                    # Check files (fasta)
            die "Fasta file ($$file_ref) has an entry!\nProgram is terminated\n";                                      # Check files (fasta)
        }                                                                                                              # Check files (fasta)
    }                                                                                                                  # Check files (fasta)
}                                                                                                                      # Check files (fasta)
## # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
sub check_read_file {                                                                                                  # Check files (read_file)
    # This subroutine check if the read files are working properly with the baiting and collecting programs            # Check files (read_file)
    #   It creates a folder to hold all the temporary files and in the end it deletes unless there was error           # Check files (read_file)
    #   Inputs:                                                                                                        # Check files (read_file)
    my ($log_ref,       # Reference to the log file handle                                                             # Check files (read_file)
        $reads_ref,     # Reference to the array holding the symlinks to the read files                                # Check files (read_file)
        $single_ref,    # Reference to the variable which is true if the reads are unpaired                            # Check files (read_file)
        $gzip_ref       # Reference to the variable that is true if the read files are compressed using gzip           # Check files (read_file)
                        ) = @_;                                                                                        # Check files (read_file)
                                                                                                                       # Check files (read_file)
    # Create folder if there is already one then it is deleted first                                                   # Check files (read_file)
    system("rm", "-r", "test_reads") if -d "test_reads";                                                               # Check files (read_file)
    mkdir "test_reads";                                                                                                # Check files (read_file)
                                                                                                                       # Check files (read_file)
    # Create test read files (get the first entry of each read file)                                                   # Check files (read_file)
    #   Replace the "N"s in the files because if the selected read has too man "N"s the baiting does not work          # Check files (read_file)
    # Create arrays to hold the input and output read files                                                            # Check files (read_file)
    my @test_reads;                                                                                                    # Check files (read_file)
    my @test_outs;                                                                                                     # Check files (read_file)
                                                                                                                       # Check files (read_file)
    # Create counter for the read files                                                                                # Check files (read_file)
    my $i;                                                                                                             # Check files (read_file)
    for my $read (@$reads_ref) {                                                                                       # Check files (read_file)
        $i++;                                                                                                          # Check files (read_file)
        # The default is that the read file is not compressed so cat can read it                                       # Check files (read_file)
        my $cat = "cat";                                                                                               # Check files (read_file)
        # Otherwise use zcat to read the compressed read file                                                          # Check files (read_file)
        $cat = "zcat" if $$gzip_ref;                                                                                   # Check files (read_file)
                                                                                                                       # Check files (read_file)
        # Get the first line of the read file to determine if it is a fasta or fastq format                            # Check files (read_file)
        my $first_line = `$cat $read | head -1`;                                                                       # Check files (read_file)
                                                                                                                       # Check files (read_file)
        # If the read file is FASTQ format then get the first entry (4 lines)                                          # Check files (read_file)
        #   Otherwise if it is FASTA then get the first two lines (that should have a sequence with sufficient length  # Check files (read_file)
        if ($first_line =~ /^\@/) {     # FASTQ                                                                        # Check files (read_file)
            # Creates an input read file                                                                               # Check files (read_file)
            `$cat $read | head -4 | sed s/N/A/g >test_reads/first_read$i\.fastq`;                                      # Check files (read_file)
        } elsif ($first_line =~ /^>/) { # FASTA                                                                        # Check files (read_file)
            # Creates an input read file                                                                               # Check files (read_file)
            `$cat $read | head -2 | sed s/N/A/g >test_reads/first_read$i\.fastq`;                                      # Check files (read_file)
        } else {                                                                                                       # Check files (read_file)
            # If the id line doesn't start either with @ or > then there is something wrong                            # Check files (read_file)
            die "The read file cannot be read or it is not in FASTA/FASTQ format!\nProgram is terminated\n";           # Check files (read_file)
        }                                                                                                              # Check files (read_file)
                                                                                                                       # Check files (read_file)
        # Add the input and output read files                                                                          # Check files (read_file)
        push @test_reads, "first_read$i\.fastq";                                                                       # Check files (read_file)
        push @test_outs, "out$i\.fastq";                                                                               # Check files (read_file)
    }                                                                                                                  # Check files (read_file)
                                                                                                                       # Check files (read_file)
    # Create test bait (the sequence of the first read in the read file)                                               # Check files (read_file)
    # The name of the test bait                                                                                        # Check files (read_file)
    my $test_bait = "test_reads/bait.fas";                                                                             # Check files (read_file)
    # Open a file handle for the bait file                                                                             # Check files (read_file)
    open my $test_bait_handle, '>', $test_bait;                                                                        # Check files (read_file)
    # Save the sequence of the first read in the read file as the bait sequence                                        # Check files (read_file)
    print {$test_bait_handle} ">first_sequence\n" . `head -2 test_reads/first_read1.fastq | tail -n 1`;                # Check files (read_file)
    # Close the file handle                                                                                            # Check files (read_file)
    close $test_bait_handle;                                                                                           # Check files (read_file)
                                                                                                                       # Check files (read_file)
    # Move to the test folder                                                                                          # Check files (read_file)
    chdir "test_reads";                                                                                                # Check files (read_file)
                                                                                                                       # Check files (read_file)
    # Create the proper relative paths to the files needed                                                             # Check files (read_file)
    $test_bait = "bait.fas";                                                                                           # Check files (read_file)
    my $test_new = "new.list";                                                                                         # Check files (read_file)
    my $test_old = "old.list";                                                                                         # Check files (read_file)
                                                                                                                       # Check files (read_file)
    # Create test collection files                                                                                     # Check files (read_file)
    system("touch", $test_old);                                                                                        # Check files (read_file)
                                                                                                                       # Check files (read_file)
    # Open a log file within the test folder for the read test                                                         # Check files (read_file)
    my $test_log = "log.log";                                                                                          # Check files (read_file)
    open my $test_handle, '>>', $test_log;                                                                             # Check files (read_file)
                                                                                                                       # Check files (read_file)
    # Run find_reads on the test files                                                                                 # Check files (read_file)
    # If the reads are not paired then test the read files individually                                                # Check files (read_file)
    if ($$single_ref) {                                                                                                # Check files (read_file)
        for (0..$#test_reads) {                                                                                        # Check files (read_file)
            # Remove files from previous test if they are present                                                      # Check files (read_file)
            unlink("new.list") if -e "new.list";                                                                       # Check files (read_file)
            unlink("positive_1.txt") if -e "positive_1.txt";                                                           # Check files (read_file)
            # Create the bait file for the current read                                                                # Check files (read_file)
	    open my $test_bait_handle, '>', $test_bait;                                                                # Check files (read_file)
	    print {$test_bait_handle} ">first_sequence\n" . `head -2 $test_reads[$_] | tail -n 1`;                     # Check files (read_file)
	    close $test_bait_handle;                                                                                   # Check files (read_file)
                                                                                                                       # Check files (read_file)
            # Create necessary local variables                                                                         # Check files (read_file)
            my @reads = ($test_reads[$_]);                                                                             # Check files (read_file)
            my @outs = ($test_outs[$_]);                                                                               # Check files (read_file)
	    my $l_tab = "";                                                                                            # Check files (read_file)
            &find_reads(\$test_handle, \$test_bait, \@reads, \$test_old, \$test_new, \@outs, \$l_tab ,"check");        # Check files (read_file)
	}                                                                                                              # Check files (read_file)
    } else {                                                                                                           # Check files (read_file)
	my $l_tab = "";                                                                                                # Check files (read_file)
	&find_reads(\$test_handle, \$test_bait, \@test_reads, \$test_old, \$test_new, \@test_outs, \$l_tab ,"check");  # Check files (read_file)
    }                                                                                                                  # Check files (read_file)
                                                                                                                       # Check files (read_file)
    # Test if the output files are correct (they should have the same line count as the input files                    # Check files (read_file)
    # Create a counter so we can refer to the proper element of the output read files array                            # Check files (read_file)
    $i = 0;                                                                                                            # Check files (read_file)
    for my $read (@test_reads) {                                                                                       # Check files (read_file)
        if ( `cat $read | wc -l` != `cat $test_outs[$i] | wc -l` ) {                                                   # Check files (read_file)
            $print = "Error: the entries in the read file must be renamed!\n";                                         # Check files (read_file)
            print {$$log_ref} $print; print $print;                                                                    # Check files (read_file)
            die "The entries in the read file must be renamed!\nProgram is terminated\n";                              # Check files (read_file)
        }                                                                                                              # Check files (read_file)
        $i++;                                                                                                          # Check files (read_file)
    }                                                                                                                  # Check files (read_file)
                                                                                                                       # Check files (read_file)
    # If the file pairs have the same length then print this to the test log file                                      # Check files (read_file)
    print {$test_handle} "\t\t\t\tThe reads could be retrieved from the same files\n";                                 # Check files (read_file)
                                                                                                                       # Check files (read_file)
    # Unless reads are unpaired we have to test if the pair of the sequences were properly identified or not           # Check files (read_file)
    unless ($$single_ref) {                                                                                            # Check files (read_file)
        # If reads are paired there should be only one read in the collection file                                     # Check files (read_file)
        if ( `cat $test_new | wc -l` != 1) {                                                                           # Check files (read_file)
            $print = "Error: the entries in the read file must be renamed (pairing is broken)!\n";                     # Check files (read_file)
            print {$$log_ref} $print; print $print;                                                                    # Check files (read_file)
            die "The entries in the read file must be renamed (pairing is broken)!\nProgram is terminated\n";          # Check files (read_file)
        }                                                                                                              # Check files (read_file)
    }                                                                                                                  # Check files (read_file)
    # Close the file handle of the test log file                                                                       # Check files (read_file)
    close $test_handle;                                                                                                # Check files (read_file)
                                                                                                                       # Check files (read_file)
    # Move back to the parent folder, same as tat the start of the subroutine                                          # Check files (read_file)
    chdir "..";                                                                                                        # Check files (read_file)
                                                                                                                       # Check files (read_file)
    # If it went okay then the folder is not necessary                                                                 # Check files (read_file)
    system("rm", "-r", "test_reads") if -d "test_reads";                                                               # Check files (read_file)
}                                                                                                                      # Check files (read_file)
########################################################################################################################
###############################################################################################################
# Find reads: baiting, merging, collecting                                                                    #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
sub find_reads{                                                                                               # Find reads
    # This subroutine collects the reads that match to the current bait or their id is already in the         # Find reads
    #   collection. The pairs of reads are collected if the reads are paired.                                 # Find reads
    #   (This is based on the naming of the read entries, so if the reads are not paired, but there are reads # Find reads
    #   with the same id then they will be collected.)                                                        # Find reads
    #   Creates a positive hit file for each read file and creates a new collection file by merging these     # Find reads
    #   then it creates the new read files for assembly with the selected reads. Finally, it compares the new # Find reads
    #   and old collection files and if there are no new reads it returns the string "No new reads"           # Find reads
    #   Returns true ("No new reads") if new read files were found                                            # Find reads
    #   Inputs:                                                                                               # Find reads
    my ($log_ref,   # Reference to the log file handle                                                        # Find reads
        $bait_ref,  # Reference to the bait file                                                              # Find reads
        $reads_ref, # Reference to the array with the input read files                                        # Find reads
        $old_ref,   # Reference to the old collection                                                         # Find reads
        $new_ref,   # Reference to the new collection                                                         # Find reads
        $readpools, # Reference to the array with the output read files                                       # Find reads
        $tab,       # Reference to the variable possibly holding some tabs                                    # Find reads
        $not_stdout # Reference to the variable that is true if out put should not be printed to the STDOUT   # Find reads
                        ) = @_;                                                                               # Find reads
                                                                                                              # Find reads
    # If $tab is undefined then put a reference to an empty string so the print commands won't throw an error # Find reads
    $tab = \"" unless $tab;#;#"                                                                               # Find reads
                                                                                                              # Find reads
    # Print the current status: Finding reads                                                                 # Find reads
    $print = "" . strftime('%H:%M:%S',localtime) . "\t$$tab\tFind matching reads\n";                          # Find reads
    print {$$log_ref} $print; print $print unless $not_stdout;                                                # Find reads
                                                                                                              # Find reads
    # The base name for the file holding the ids of the positive reads                                        # Find reads
    my $mira_temp = "positive";                                                                               # Find reads
                                                                                                              # Find reads
    # For each read file do mirabait and store the positive list files in an array                            # Find reads
    my @hit_files;                                                                                            # Find reads
    my $i;                                                                                                    # Find reads
    for my $read (@$reads_ref) {                                                                              # Find reads
        $i++;                                                                                                 # Find reads
        # Check if there is already a non-empty mira output file or not                                       # Find reads
        if (not -e "$mira_temp\_$i\.txt" || -z "$mira_temp\_$i\.txt") {                                       # Find reads
            # Run mirabait                                                                                    # Find reads
            #   mirabait -t txt bait_file read_file output_prefix                                             # Find reads
            system("$bait_cmd -t txt $$bait_ref $read $mira_temp\_$i 2>&1 >>mirabait.log");                   # Find reads
	    # Test if baiting program is fuctioning correctly                                                 # Find reads
	    unless ($? == 0) {                                                                                # Find reads
		# bait has crashed                                                                            # Find reads
		$print = "Error: Baiting failed!\n";                                                          # Find reads
		$print .= "\tThe command that crashed is:\n";                                                 # Find reads
		$print .= "\t\t' $bait_cmd -t txt $$bait_ref $read $mira_temp\_$i'\n";                        # Find reads
		$print .= "\tIf you are using mirabait, try installing another version (e.g. 4.0).\n";        # Find reads
		$print .= "\t\tMira 4.0.2 has a bug and will crash when using the above command\n";           # Find reads
		unless ($not_stdout) {                                                                        # Find reads
		    print {$$log_ref} $print; print $print;                                                   # Find reads
		} else {                                                                                      # Find reads
		    print {$$log_ref} $print;                                                                 # Find reads
		    die $print;                                                                               # Find reads
		}                                                                                             # Find reads
	    }                                                                                                 # Find reads
            $print = strftime('%H:%M:%S',localtime) . "\t$$tab\t\tDone baiting read file #$i\n";              # Find reads
            print {$$log_ref} $print; print $print unless $not_stdout;                                        # Find reads
        }                                                                                                     # Find reads
                                                                                                              # Find reads
        # If there is no positive file then create an empty file, so the program doesn't crash                # Find reads
        unless (-e "$mira_temp\_$i\.txt") {                                                                   # Find reads
            $print = "Error: mirabiat failed, no $mira_temp\_$i\.txt\n";                                      # Find reads
            print {$$log_ref} $print; print $print unless $not_stdout;                                        # Find reads
            open my $empty, '>>', "$mira_temp\_$i\.txt";                                                      # Find reads
            close $empty;                                                                                     # Find reads
            $print = "\t\t$$tab\tCreated empty file $mira_temp\_$i\.txt\n";                                   # Find reads
            print {$$log_ref} $print; print $print unless $not_stdout;                                        # Find reads
        }                                                                                                     # Find reads
        # Store away the positive list files name in the array                                                # Find reads
        push @hit_files, "$mira_temp\_$i\.txt";                                                               # Find reads
    }                                                                                                         # Find reads
                                                                                                              # Find reads
    # Merge old collection and the new mirabait results THEN create new collection                            # Find reads
    # Skip if new collection already exists                                                                   # Find reads
    if (-e $$new_ref) {                                                                                       # Find reads
        $print = "" . strftime('%H:%M:%S',localtime) . "\t$$tab\tThe found merged file is used\n";            # Find reads
        print {$$log_ref} $print; print $print unless $not_stdout;                                            # Find reads
    } else {                                                                                                  # Find reads
        $print = "" . strftime('%H:%M:%S',localtime) . "\t$$tab\tMerge new hits with old collection\n";       # Find reads
        print {$$log_ref} $print; print $print unless $not_stdout;                                            # Find reads
        &merge_collections(\@hit_files, $old_ref, $new_ref);                                                  # Find reads
    }                                                                                                         # Find reads
                                                                                                              # Find reads
    # IF there are no new reads THEN stop                                                                     # Find reads
    #   ELSE continue                                                                                         # Find reads
    # Returns "No new reads" if there were no new reads found                                                 # Find reads
    my $no_new_reads = &compare_collections($log_ref, $old_ref, $new_ref, $tab, $not_stdout);                 # Find reads
    if ($no_new_reads) {                                                                                      # Find reads
        return $no_new_reads;                                                                                 # Find reads
    }                                                                                                         # Find reads
                                                                                                              # Find reads
    # Create new readpools (output read files)                                                                # Find reads
    $print = "" . strftime('%H:%M:%S',localtime) . "\t$$tab\tCreate new read files\n";                        # Find reads
    print {$$log_ref} $print; print $print unless $not_stdout;                                                # Find reads
                                                                                                              # Find reads
    # Run collect reads for each non-existing readpool file                                                   # Find reads
    #   Use a counter to refer to the right read file from the array                                          # Find reads
    $i = 0;                                                                                                   # Find reads
    for my $readpool (@$readpools) {                                                                          # Find reads
        &collect_reads(\$reads_ref->[$i], $new_ref, \$readpool) unless -e $readpool;                          # Find reads
        $i++;                                                                                                 # Find reads
    }                                                                                                         # Find reads
                                                                                                              # Find reads
    return;                                                                                                   # Find reads
}                                                                                                             # Find reads
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
sub merge_collections{                                                                                        # Find reads (merge collections)
    # This subroutine takes a few list files and merges them into a single file                               # Find reads (merge collections)
    #   Inputs:                                                                                               # Find reads (merge collections)
    my ($hits_ref,  # Reference to an array of the positive hit lists                                         # Find reads (merge collections)
    $old_ref,       # Reference to the old collection file                                                    # Find reads (merge collections)
    $new_ref        # Reference to the new read file                                                          # Find reads (merge collections)
                    ) = @_;                                                                                   # Find reads (merge collections)
    # Merge the files and remove duplicate lines                                                              # Find reads (merge collections)
    #   Remove possible pair markers from the end of the ids                                                  # Find reads (merge collections)
    `cat @$hits_ref $$old_ref | perl -ne 's/\\/1\$//; s/\\/2\$//; print' | sort | uniq >$$new_ref`;           # Find reads (merge collections)
                                                                                                              # Find reads (merge collections)
}                                                                                                             # Find reads (merge collections)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
sub collect_reads{                                                                                            # Find reads (collect reads)
    # This subroutine collects the reads that have their ids in the collection file and writes them in a read # Find reads (collect reads)
    #   file for assembly for each input read file.                                                           # Find reads (collect reads)
    #   Inputs:                                                                                               # Find reads (collect reads)
    my ($read_ref,      # Reference to the input read file                                                    # Find reads (collect reads)
        $list_ref,      # Reference to the collection file                                                    # Find reads (collect reads)
        $readpool_ref   # Reference to the output read file (for the assembly)                                # Find reads (collect reads)
                        ) = @_;                                                                               # Find reads (collect reads)
                                                                                                              # Find reads (collect reads)
    # Collect the reads from the input read file into the output read files                                   # Find reads (collect reads)
    #   Use the specified collecting program to collect the reads                                             # Find reads (collect reads)
    `$collect_cmd $$read_ref $$list_ref >$$readpool_ref`;                                                     # Find reads (collect reads)
}                                                                                                             # Find reads (collect reads)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
sub compare_collections{                                                                                      # Find reads (check collections)
    # This subroutine compares two collection files based on number of non-empty lines.                       # Find reads (check collections)
    #   If the new collection has more reads then it returns FALSE else it returns "No new reads"             # Find reads (check collections)
    #   Inputs:                                                                                               # Find reads (check collections)
    my ($log_ref,   # Reference to the log file handle                                                        # Find reads (check collections)
        $old_ref,   # Reference to the old collection                                                         # Find reads (check collections)
        $new_ref,   # Reference to the new collection                                                         # Find reads (check collections)
        $tab,       # Reference to the variable possibly holding some tabs                                    # Find reads (check collections)
        $not_stdout # Reference to the variable that is true if out put should not be printed to the STDOUT   # Find reads (check collections)
                    ) = @_;                                                                                   # Find reads (check collections)
                                                                                                              # Find reads (check collections)
    # Print the current status: checking read count                                                           # Find reads (check collections)
    $print = "" . strftime('%H:%M:%S',localtime) . "\t$$tab\tChecking read count\n";                          # Find reads (check collections)
    print {$$log_ref} $print; print $print unless $not_stdout;                                                # Find reads (check collections)
                                                                                                              # Find reads (check collections)
    # Create a variable for each file, that will hold the line number                                         # Find reads (check collections)
    my ($oldC, $newC);                                                                                        # Find reads (check collections)
                                                                                                              # Find reads (check collections)
    # IF old collection does NOT EXIST make the previous read count 0 and mention the missing file in the log # Find reads (check collections)
    #   ELSE get the read count                                                                               # Find reads (check collections)
    unless (-e $$old_ref) {                                                                                   # Find reads (check collections)
        $oldC = 0;                                                                                            # Find reads (check collections)
        $print = "\t\t\t$$tab\tOld collection file ($$old_ref) does not exist\n";                             # Find reads (check collections)
        print {$$log_ref} $print; print $print unless $not_stdout;                                            # Find reads (check collections)
    } else {                                                                                                  # Find reads (check collections)
        # Count the number of non-empty lines in the collection                                               # Find reads (check collections)
        $oldC = `cat $$old_ref | sed '/^\\s*\$/d' | wc -l`;                                                   # Find reads (check collections)
    }                                                                                                         # Find reads (check collections)
                                                                                                              # Find reads (check collections)
    # IF new collection does NOT EXIST then there are 0 ids in it                                             # Find reads (check collections)
    #   ELSE get the read count from the file                                                                 # Find reads (check collections)
    unless (-e $$new_ref) {                                                                                   # Find reads (check collections)
        $print = "ERROR: New collection file ($$new_ref) does not exist\n";                                   # Find reads (check collections)
        print {$$log_ref} $print; print $print unless $not_stdout;                                            # Find reads (check collections)
        $newC = 0;                                                                                            # Find reads (check collections)
    } else {                                                                                                  # Find reads (check collections)
        # Count the number of non-empty lines in the collection                                               # Find reads (check collections)
        $newC = `cat $$new_ref | sed '/^\\s*\$/d' | wc -l`;                                                   # Find reads (check collections)
    }                                                                                                         # Find reads (check collections)
                                                                                                              # Find reads (check collections)
    # Stop if no reads are found                                                                              # Find reads (check collections)
    unless ($newC - 0) {                                                                                      # Find reads (check collections)
        $print = "" . strftime('%H:%M:%S',localtime) . "\t\t$$tab" . "No reads were found\n";                 # Find reads (check collections)
        print {$$log_ref} $print; print $print unless $not_stdout;                                            # Find reads (check collections)
        # Differentiate between no reads and no new reads                                                     # Find reads (check collections)
        return "No reads";                        # TRUE                                                      # Find reads (check collections)
    }                                                                                                         # Find reads (check collections)
    # RETURN TRUE UNLESS the second file has MORE non-empty lines than the first ELSE FALSE                   # Find reads (check collections)
    unless ($oldC < $newC) {                                                                                  # Find reads (check collections)
        # Because the read count did not increase there cannot be new reads                                   # Find reads (check collections)
        $print = "" . strftime('%H:%M:%S',localtime) . "\t\t$$tab\tNumber of new reads:   0\n";               # Find reads (check collections)
        print {$$log_ref} $print; print $print unless $not_stdout;                                            # Find reads (check collections)
        $print = strftime('%H:%M:%S',localtime) . "\t\t$$tab\tNumber of total reads: " . ($oldC - 0) . "\n";  # Find reads (check collections)
        print {$$log_ref} $print; print $print unless $not_stdout;                                            # Find reads (check collections)
        return "No new reads"; 	                  # TRUE                                                      # Find reads (check collections)
    } else {                                                                                                  # Find reads (check collections)
        # Print the read counts unless STDOUT is suppressed                                                   # Find reads (check collections)
        $print = "" . strftime('%H:%M:%S',localtime);                                                         # Find reads (check collections)
        print {$$log_ref} $print; print $print unless $not_stdout;                                            # Find reads (check collections)
        $print = "\t\t$$tab\tNumber of new reads:   " . ($newC - $oldC) . "\n";                               # Find reads (check collections)
        print {$$log_ref} $print; print $print unless $not_stdout;                                            # Find reads (check collections)
        $print = "" . strftime('%H:%M:%S',localtime);                                                         # Find reads (check collections)
        print {$$log_ref} $print; print $print unless $not_stdout;                                            # Find reads (check collections)
        $print = "\t\t$$tab\tNumber of total reads: " . ($newC - 0) . "\n";                                   # Find reads (check collections)
        print {$$log_ref} $print; print $print unless $not_stdout;                                            # Find reads (check collections)
        return;                 # FALSE                                                                       # Find reads (check collections)
    }                                                                                                         # Find reads (check collections)
}                                                                                                             # Find reads (check collections)
###############################################################################################################
#######################################################################################################################
# Do assembly: graph generation, actual assembly                                                                      #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
sub assemble{                                                                                                         # Assemble
    # Input: reads[array(filename strings)], parameters[hash(key-value)], name_of_the_assembly_file[string]           # Assemble
    #   whatever the assembler outputs as final assembly it is copied to $assembly_output                             # Assemble
    #   Inputs:                                                                                                       # Assemble
    my ($log_ref,           # Reference to the log file handle                                                        # Assemble
        $readpool_ref,      # Reference to the array that holds the read file names for the assembly                  # Assemble
        $single_ref,        # Reference to the variable that is true if the reads are unpaired                        # Assemble
        $parameters_ref,    # Reference to an array that holds all the arguments that needed for the assembly         # Assemble
                            #   [0] name of the assembler program                                                     # Assemble
                            #   [1] additional parameters for the graph/hash generation                               # Assemble
                            #   [2] additional parameters for the assembly generation                                 # Assemble
        $out_file_ref,      # Reference to the file name for the assembly                                             # Assemble
        $tab,               # Reference to the variable possibly holding some tabs                                    # Assemble
        $format_ref,        # Reference to the variable that holds the read file format                               # Assemble
                            ) = @_;                                                                                   # Assemble
                                                                                                                      # Assemble
    # If $tab was not defined then use a reference to an empty string                                                 # Assemble
    $tab = \"" unless $tab;#"                                                                                         # Assemble
    # Print the current state of the program                                                                          # Assemble
    $print = "" . strftime('%H:%M:%S',localtime) . "\t$$tab\tAssemble collected reads\n";                             # Assemble
    print {$$log_ref} $print; print $print;                                                                           # Assemble
                                                                                                                      # Assemble
    # If the specified assembly file already exists skip assembly and print that it was skipped                       # Assemble
    #   Otherwise check which assembler subroutine should be used                                                     # Assemble
    if (-e $$out_file_ref) {                                                                                          # Assemble
        $print = "\t\t$$tab\tUsing assembly file from previous run\n";                                                # Assemble
        print {$$log_ref} $print; print $print;                                                                       # Assemble
    } else {                                                                                                          # Assemble
        # If the first element of the parameter array is edena then use the edena assembler subroutine                # Assemble
        if ($$parameters_ref[0] eq "edena") {                                                                         # Assemble
            &edena_assemble($log_ref, $readpool_ref, $single_ref, $parameters_ref, $out_file_ref, $tab);              # Assemble
            # If the first element of the parameter array is alternative then use the alternative subroutine          # Assemble
        } elsif ($$parameters_ref[0] eq "alternative") {                                                              # Assemble
            &alternative_assemble($log_ref, $readpool_ref, $single_ref, $parameters_ref,                              # Assemble
                                  $out_file_ref, $tab, $format_ref);                                                  # Assemble
            # If the first element of the parameter array is velvet then use the velvet assembler subroutine          # Assemble
        } elsif ($$parameters_ref[0] eq "velvet") {                                                                   # Assemble
            &velvet_assemble($log_ref, $readpool_ref, $single_ref, $parameters_ref, $out_file_ref, $tab, $format_ref);# Assemble
        } elsif ($$parameters_ref[0] eq "external") {                                                                 # Assemble
            my $arg_single = "paired";                                                                                # Assemble
            $arg_single = "single" if $$single_ref;                                                                   # Assemble
            my $extra1 = "-";                                                                                         # Assemble
            my $extra2 = "-";                                                                                         # Assemble
            $extra1 = $$parameters_ref[1] if $$parameters_ref[1];                                                     # Assemble
            $extra2 = $$parameters_ref[2] if $$parameters_ref[2];                                                     # Assemble
            `perl ../../assembler.pl "@$readpool_ref" $arg_single "$extra1" "$extra2" $$out_file_ref $$format_ref`;   # Assemble
        }                                                                                                             # Assemble
        #                                                                                                             # Assemble
        # If there is no assembly file print an error                                                                 # Assemble
        #   Create an empty file instead, maybe the other threads are still running                                   # Assemble
        unless (-e $$out_file_ref) {                                                                                  # Assemble
            $print = "Error: Assembly failed, output file ($$out_file_ref) is missing.\n";                            # Assemble
            print {$$log_ref} $print; print $print;                                                                   # Assemble
            #                                                                                                         # Assemble
            # Create an empty file instead of the output file to avoid crashing the program                           # Assemble
            open my $empty , '>>', $$out_file_ref;                                                                    # Assemble
            close $empty;                                                                                             # Assemble
            $print = "\t\t$$tab\tCreated an empty file ($$out_file_ref)\n";                                           # Assemble
            print {$$log_ref} $print; print $print;                                                                   # Assemble
        }                                                                                                             # Assemble
    }                                                                                                                 # Assemble
}                                                                                                                     # Assemble
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
sub edena_assemble{                                                                                                   # Assemble (Edena)
    # This subroutine uses the edena assembler to create the assembly. The assembly file at the end is copied to the  # Assemble (Edena)
    #   specified file ($out_file_ref)                                                                                # Assemble (Edena)
    #   Inputs:                                                                                                       # Assemble (Edena)
    my ($log_ref,           # Reference to the log file handle                                                        # Assemble (Edena)
    $readpool_ref,          # Reference to the array that holds the read file names for the assembly                  # Assemble (Edena)
    $single_ref,            # Reference to the variable that is true if the reads are unpaired                        # Assemble (Edena)
    $parameters_ref,        # Reference to an array that holds all the arguments that needed for the assembly         # Assemble (Edena)
                            #   [0] name of the assembler program                                                     # Assemble (Edena)
                            #   [1] additional parameters for the graph/hash generation                               # Assemble (Edena)
                            #   [2] additional parameters for the assembly generation                                 # Assemble (Edena)
    $out_file_ref,          # Reference to the file name for the assembly                                             # Assemble (Edena)
    $tab                    # Reference to the variable possibly holding some tabs                                    # Assemble (Edena)
                            ) = @_;                                                                                   # Assemble (Edena)
                                                                                                                      # Assemble (Edena)
    # If $tab was not defined then use a reference to an empty string                                                 # Assemble (Edena)
    $tab = \"" unless $tab;#"                                                                                         # Assemble (Edena)
                                                                                                                      # Assemble (Edena)
    # pairing specific input for edena                                                                                # Assemble (Edena)
    my $pairing;                                                                                                      # Assemble (Edena)
    if ($$single_ref) {                                                                                               # Assemble (Edena)
        $pairing = "-singleEnd";                                                                                      # Assemble (Edena)
    } else {                                                                                                          # Assemble (Edena)
        $pairing = "-paired";                                                                                         # Assemble (Edena)
    }                                                                                                                 # Assemble (Edena)
                                                                                                                      # Assemble (Edena)
    #  Print the status of the program: Overlap graph generation                                                      # Assemble (Edena)
    $print = strftime('%H:%M:%S',localtime) . "\t\t$$tab\tRunning overlap analysis using Edena assembler\n";          # Assemble (Edena)
    print {$$log_ref} $print; print $print;                                                                           # Assemble (Edena)
                                                                                                                      # Assemble (Edena)
    # @arg1 will hold the parameters the were specified for the graph/hash generation                                 # Assemble (Edena)
    my @arg1;                                                                                                         # Assemble (Edena)
    if ($$parameters_ref[1]) {                                                                                        # Assemble (Edena)
        @arg1 = split / /, $$parameters_ref[1];                                                                       # Assemble (Edena)
    }                                                                                                                 # Assemble (Edena)
                                                                                                                      # Assemble (Edena)
    # Command line command for graph/hash generation                                                                  # Assemble (Edena)
    `$edena_cmd $pairing @$readpool_ref @arg1 -prefix edena 2>&1 >>edena_graph.log`;                                  # Assemble (Edena)
                                                                                                                      # Assemble (Edena)
                                                                                                                      # Assemble (Edena)
    # If there is no overlap file then the assembly has failed, so return                                             # Assemble (Edena)
    unless (-e "edena.ovl") {                                                                                         # Assemble (Edena)
        $print = "Error: Assembly failed, overlap file (edena.ovl) is missing\n";                                     # Assemble (Edena)
        print {$$log_ref} $print; print $print;                                                                       # Assemble (Edena)
        return;                                                                                                       # Assemble (Edena)
    }                                                                                                                 # Assemble (Edena)
                                                                                                                      # Assemble (Edena)
    # Print the status of the program: Actual assembly                                                                # Assemble (Edena)
    $print = strftime('%H:%M:%S',localtime) . "\t\t$$tab\tRunning assembly using Edena assembler\n";                  # Assemble (Edena)
    print {$$log_ref} $print; print $print;                                                                           # Assemble (Edena)
                                                                                                                      # Assemble (Edena)
    # @arg2 will hold the parameters the were specified for the assembly generation                                   # Assemble (Edena)
    my @arg2;                                                                                                         # Assemble (Edena)
    if ($$parameters_ref[2]) {                                                                                        # Assemble (Edena)
        @arg2 = split / /, $$parameters_ref[2];                                                                       # Assemble (Edena)
    }                                                                                                                 # Assemble (Edena)
                                                                                                                      # Assemble (Edena)
    # Command line command for assembly                                                                               # Assemble (Edena)
    `$edena_cmd -e edena.ovl @arg2 -prefix edena 2>&1 >>edena_asmbl.log`;                                             # Assemble (Edena)
                                                                                                                      # Assemble (Edena)
    # Print if there is no assembly file generated                                                                    # Assemble (Edena)
    unless (-e "edena_contigs.fasta") {                                                                               # Assemble (Edena)
        $print = "Error: Assembly failed, contig file (edena_contigs.fasta) is missing\n";                            # Assemble (Edena)
        print {$$log_ref} $print; print $print;                                                                       # Assemble (Edena)
        return;                                                                                                       # Assemble (Edena)
    } else {                                                                                                          # Assemble (Edena)
        # Copy the assembly file to the specified file                                                                # Assemble (Edena)
        system("cp", "edena_contigs.fasta", $$out_file_ref);                                                          # Assemble (Edena)
    }                                                                                                                 # Assemble (Edena)
}                                                                                                                     # Assemble (Edena)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
sub velvet_assemble{                                                                                                  # Assemble (Velvet)
    # This subroutine uses the velvet assembler to create the assembly. The assembly file at the end is copied to the # Assemble (Velvet)
    #   specified file ($out_file_ref)                                                                                # Assemble (Velvet)
    #   Inputs:                                                                                                       # Assemble (Velvet)
    my ($log_ref,           # Reference to the log file handle                                                        # Assemble (Velvet)
        $readpool_ref,      # Reference to the array that holds the read file names for the assembly                  # Assemble (Velvet)
        $single_ref,        # Reference to the variable that is true if the reads are unpaired                        # Assemble (Velvet)
        $parameters_ref,    # Reference to an array that holds all the arguments that needed for the assembly         # Assemble (Velvet)
                            #   [0] name of the assembler program                                                     # Assemble (Velvet)
                            #   [1] additional parameters for the graph/hash generation                               # Assemble (Velvet)
                            #   [2] additional parameters for the assembly generation                                 # Assemble (Velvet)
        $out_file_ref,      # Reference to the file name for the assembly                                             # Assemble (Velvet)
        $tab,               # Reference to the variable possibly holding some tabs                                    # Assemble (Velvet)
        $format_ref,        # Reference to the variable that holds the read file format                               # Assemble (Velvet)
                            ) = @_;                                                                                   # Assemble (Velvet)
                                                                                                                      # Assemble (Velvet)
    # If $tab was not defined then use a reference to an empty string                                                 # Assemble (Velvet)
    $tab = \"" unless $tab;#"                                                                                         # Assemble (Velvet)
                                                                                                                      # Assemble (Velvet)
    # @arg1 will hold the parameters the were specified for the graph/hash generation                                 # Assemble (Velvet)
    my @arg1;                                                                                                         # Assemble (Velvet)
    if ($$parameters_ref[1]) {                                                                                        # Assemble (Velvet)
        @arg1 = split / /, $$parameters_ref[1];                                                                       # Assemble (Velvet)
    }                                                                                                                 # Assemble (Velvet)
                                                                                                                      # Assemble (Velvet)
    # @arg2 will hold the parameters the were specified for the assembly generation                                   # Assemble (Velvet)
    my @arg2;                                                                                                         # Assemble (Velvet)
    if ($$parameters_ref[2]) {                                                                                        # Assemble (Velvet)
        @arg2 = split / /, $$parameters_ref[2];                                                                       # Assemble (Velvet)
    }                                                                                                                 # Assemble (Velvet)
                                                                                                                      # Assemble (Velvet)
    # Name of the assembler that should be printed                                                                    # Assemble (Velvet)
    my $name = "Velvet"; # name of the assembler (using Velvet as an example)                                         # Assemble (Velvet) !!!!
    # Graph generation using single-end mode (@arg1 holds the argument specified at invocation)                       # Assemble (Velvet)
    my $single_first = "$velveth_cmd velvet 31 -short -$$format_ref @$readpool_ref @arg1";                            # Assemble (Velvet) !!!!
    # Assembly using single-end mode (@arg2 holds the argument specified at invocation)                               # Assemble (Velvet)
    my $single_second = "$velvetg_cmd velvet @arg2"; # The command line command for actual assembly                   # Assemble (Velvet) !!!!
    # Graph generation using paired-end mode (@arg1 holds the argument specified at invocation)                       # Assemble (Velvet)
    # command is : velveth velvet 31 -short -fastq readpool1.fq readpool2.fq                                          # Assemble (Velvet) !!!!
    my $paired_first = "$velveth_cmd velvet 31 -short -$$format_ref @$readpool_ref @arg1";                            # Assemble (Velvet) !!!!
    # Assembly using paired-end mode (@arg2 holds the argument specified at invocation)                               # Assemble (Velvet)
    my $paired_second = "$velvetg_cmd velvet @arg2"; # The command line command for actual assembly                   # Assemble (Velvet) !!!!
    my $result = "velvet/contigs.fa";                                                                                 # Assemble (Velvet) !!!!
                                                                                                                      # Assemble (Velvet)
    # Need to run a different commands for single and for paired reads                                                # Assemble (Velvet)
    if ($$single_ref) {                                                                                               # Assemble (Velvet)
        # Single mode                                                                                                 # Assemble (Velvet)
        # Overlap graph generation                                                                                    # Assemble (Velvet)
        $print = strftime('%H:%M:%S',localtime) . "\t\t$$tab\tRunning overlap analysis using $name assembler\n";      # Assemble (Velvet)
        print {$$log_ref} $print; print $print;                                                                       # Assemble (Velvet)
                                                                                                                      # Assemble (Velvet)
        # Command line command for graph/hash generation                                                              # Assemble (Velvet)
        `$single_first`;                                                                                              # Assemble (Velvet)
                                                                                                                      # Assemble (Velvet)
        # Actual assembly                                                                                             # Assemble (Velvet)
        $print = "" . strftime('%H:%M:%S',localtime) . "\t\t$$tab\tRunning assembly using $name assembler\n";         # Assemble (Velvet)
        print {$$log_ref} $print; print $print;                                                                       # Assemble (Velvet)
                                                                                                                      # Assemble (Velvet)
        # Command line command for assembly generation                                                                # Assemble (Velvet)
        `$single_second`;                                                                                             # Assemble (Velvet)
    } else {                                                                                                          # Assemble (Velvet)
        # Paired-end mode                                                                                             # Assemble (Velvet)
        # Overlap graph generation                                                                                    # Assemble (Velvet)
        $print = strftime('%H:%M:%S',localtime) . "\t\t$$tab\tRunning overlap analysis using $name assembler\n";      # Assemble (Velvet)
        print {$$log_ref} $print; print $print;                                                                       # Assemble (Velvet)
                                                                                                                      # Assemble (Velvet)
        # Command line command for graph/hash generation                                                              # Assemble (Velvet)
        `$paired_first`;                                                                                              # Assemble (Velvet)
                                                                                                                      # Assemble (Velvet)
        # Actual assembly                                                                                             # Assemble (Velvet)
        $print = "" . strftime('%H:%M:%S',localtime) . "\t\t$$tab\tRunning assembly using $name assembler\n";         # Assemble (Velvet)
        print {$$log_ref} $print; print $print;                                                                       # Assemble (Velvet)
                                                                                                                      # Assemble (Velvet)
        # Command line command for assembly generation                                                                # Assemble (Velvet)
        `$paired_second`;                                                                                             # Assemble (Velvet)
    }                                                                                                                 # Assemble (Velvet)
                                                                                                                      # Assemble (Velvet)
    # Copy the assembly (output) file to the specified location                                                       # Assemble (Velvet)
    # Specify the output file of the assembler                                                                        # Assemble (Velvet)
    #   use relative path by which it can be found from where the assembly commands are run                           # Assemble (Velvet)
    my $result = "velvet/contigs.fa";                                                                                 # Assemble (Velvet) !!!!
                                                                                                                      # Assemble (Velvet)
    system("cp", $result, $$out_file_ref) if -e $result;                                                              # Assemble (Velvet)
}                                                                                                                     # Assemble (Velvet)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# A subroutine to specify a new program as assembler                                                                  # README
#   Check velvet_assemble to see how it may look like                                                                 # README
#   Don't forget to:	To specify --assembler alternative at invocation                                              # README
#   Fill in the values in the lines below marked by exclamation marks at the line ends (follow instruction)           # README
#       1) Specify the name of the assembler                                                                          # README
#       2) Write the command for graph generation (both in paired-end and single-end mode)                            # README
#       3) Write the command for assembly generation (both in paired-end and single-end mode)                         # README
#       4) Specify the output file of the assembler                                                                   # README
sub alternative_assemble{                                                                                             # Assemble (Alternative)
    # This subroutine uses the user defined assembler to do the assembly, the marked lines should be changed to make  # Assemble (Alternative)
    #   it operational                                                                                                # Assemble (Alternative)
    #   Inputs:                                                                                                       # Assemble (Alternative)
    my ($log_ref,           # Reference to the log file handle                                                        # Assemble (Alternative)
        $readpool_ref,      # Reference to the array that holds the read file names for the assembly                  # Assemble (Alternative)
        $single_ref,        # Reference to the variable that is true if the reads are unpaired                        # Assemble (Alternative)
        $parameters_ref,    # Reference to an array that holds all the arguments that needed for the assembly         # Assemble (Alternative)
                            #   [0] name of the assembler program                                                     # Assemble (Alternative)
                            #   [1] additional parameters for the graph/hash generation                               # Assemble (Alternative)
                            #   [2] additional parameters for the assembly generation                                 # Assemble (Alternative)
        $out_file_ref,      # Reference to the file name for the assembly                                             # Assemble (Alternative)
        $tab,               # Reference to the variable possibly holding some tabs                                    # Assemble (Alternative)
        $format_ref,        # Reference to the variable that holds the read file format                               # Assemble (Alternative)
                            ) = @_;                                                                                   # Assemble (Alternative)
                                                                                                                      # Assemble (Alternative)
    # If $tab was not defined then use a reference to an empty string                                                 # Assemble (Alternative)
    $tab = \"" unless $tab;#"                                                                                         # Assemble (Alternative)
                                                                                                                      # Assemble (Alternative)
    # @arg1 will hold the parameters the were specified for the graph/hash generation                                 # Assemble (Alternative)
    my @arg1;                                                                                                         # Assemble (Alternative)
    if ($$parameters_ref[1]) {                                                                                        # Assemble (Alternative)
        @arg1 = split / /, $$parameters_ref[1];                                                                       # Assemble (Alternative)
    }                                                                                                                 # Assemble (Alternative)
    # @arg2 will hold the parameters the were specified for the assembly generation                                   # Assemble (Alternative)
    my @arg2;                                                                                                         # Assemble (Alternative)
    if ($$parameters_ref[2]) {                                                                                        # Assemble (Alternative)
        @arg2 = split / /, $$parameters_ref[2];                                                                       # Assemble (Alternative)
    }                                                                                                                 # Assemble (Alternative)
                                                                                                                      # Assemble (Alternative)
    # $$format_ref is the format of the read files (either "fastq" or "fasta")                                        # Assemble (Alternative)
                                                                                                                      # Assemble (Alternative)
    #          vvvvvvvvvvvv                                                                                           # Assemble (Alternative)
    # Name of the assembler that should be printed                                                                    # Assemble (Alternative)
    my $name = ""; # name of the assembler, this will be print to the log and STDOUT    # Assemble (Alternative) !!!!!!!!!!
    # Graph generation using single-end mode (@arg1 holds the argument specified at invocation)                       # Assemble (Alternative)
    my $single_first = "$alt_1_cmd @arg1 @$readpool_ref"; # Assemble (Alternative) !!!!!!!!!
    # Assembly using single-end mode (@arg2 holds the argument specified at invocation)                               # Assemble (Alternative)
    my $single_second = "$alt_2_cmd @arg2"; # The command line command for actual assembly  # Assemble (Alternative) !!!!!!!!!
    # Graph generation using paired-end mode (@arg1 holds the argument specified at invocation)                       # Assemble (Alternative)
    my $paired_first = "$alt_1_cmd @arg1 @$readpool_ref"; # Assemble (Alternative) !!!!!!!!!
    # Assembly using paired-end mode (@arg2 holds the argument specified at invocation)                               # Assemble (Alternative)
    my $paired_second = "$alt_2_cmd @arg2"; # The command line command for actual assembly  # Assemble (Alternative) !!!!!!!!!
    my $result = "";  # relative path of the output file generated by the assembler            # Assemble (Alternative) !!!!!!!!!!!!!!!!!!!!!
    #          ^^^^^^^^^^^^                                                                                           # Assemble (Alternative)
                                                                                                                      # Assemble (Alternative)
    # Need to run a different commands for single and for paired reads                                                # Assemble (Alternative)
    if ($$single_ref) {                                                                                               # Assemble (Alternative)
        # Single mode                                                                                                 # Assemble (Alternative)
        # Overlap graph generation                                                                                    # Assemble (Alternative)
        $print = strftime('%H:%M:%S',localtime) . "\t\t\tRunning overlap analysis using $name assembler\n";           # Assemble (Alternative)
        print {$$log_ref} $print; print $print;                                                                       # Assemble (Alternative)
                                                                                                                      # Assemble (Alternative)
        `$single_first`;                                                                                              # Assemble (Alternative)
                                                                                                                      # Assemble (Alternative)
        # Actual assembly                                                                                             # Assemble (Alternative)
        $print = strftime('%H:%M:%S',localtime) . "\t\t\tRunning assembly using $name assembler\n";                   # Assemble (Alternative)
        print {$$log_ref} $print; print $print;                                                                       # Assemble (Alternative)
                                                                                                                      # Assemble (Alternative)
        `$single_second`;                                                                                             # Assemble (Alternative)
                                                                                                                      # Assemble (Alternative)
    } else {                                                                                                          # Assemble (Alternative)
        # Paired-end mode                                                                                             # Assemble (Alternative)
        # Overlap graph generation                                                                                    # Assemble (Alternative)
        $print = strftime('%H:%M:%S',localtime) . "\t\t\tRunning overlap analysis using $name assembler\n";           # Assemble (Alternative)
        print {$$log_ref} $print; print $print;                                                                       # Assemble (Alternative)
                                                                                                                      # Assemble (Alternative)
        `$paired_first`;                                                                                              # Assemble (Alternative)
                                                                                                                      # Assemble (Alternative)
        # Actual assembly                                                                                             # Assemble (Alternative)
        $print = strftime('%H:%M:%S',localtime) . "\t\t\tRunning assembly using $name assembler\n";                   # Assemble (Alternative)
        print {$$log_ref} $print; print $print;                                                                       # Assemble (Alternative)
                                                                                                                      # Assemble (Alternative)
        `$paired_second`;                                                                                             # Assemble (Alternative) 
                                                                                                                      # Assemble (Alternative)
    }                                                                                                                 # Assemble (Alternative)
    # Copy the assembly (output) file to the specified location                                                       # Assemble (Alternative)
    # Specify the output file of the assembler                                                                        # Assemble (Alternative)
    #   use relative path by which it can be found from where the assembly commands are run                           # Assemble (Alternative)
                                                                                                                      # Assemble (Alternative)
    system("cp", $result, $$out_file_ref) if -e $result;                                                              # Assemble (Alternative)
}                                                                                                                     # Assemble (Alternative)
#######################################################################################################################
###############################################################################################################
# Test if run is completed                                                                                    #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
sub test_completion{                                                                                          # Test completion
    ## Check if another iteration is needed or not                                                            # Test completion
    # IF it is complete THEN RETURN true                                                                      # Test completion
    #   These are the different criteria and their return/true values                                         # Test completion
    #   Criteria:    1) No or empty bait                    =>"No assembly"                                   # Test completion
    #                2) Bait has no new information         =>"Same bait"                                     # Test completion
    #                3) Using exonerate to find homologous region                                             # Test completion
    #                    matches to a single contig         =>"Found"                                         # Test completion
    #                    matches to multiple contigs        =>"Multiple contigs"                              # Test completion
    #                    no improvement in coverage         =>"No improvement"                                # Test completion
    #                4) The total contig size               =>"Reached: total length is at least $min bp"     # Test completion
    #                5) The length of the longest contig    =>"Reached: longest contig is at least $min bp"   # Test completion
    #                6) The N50 of the assembly             =>"Reached: N50 is at least $min bp"              # Test completion
    #       The first two are exclusive and also excludes the possibility of the next criteria being true     # Test completion
    #       criteria 3-6 are not exclusive so if multiple criteria are used then the first met criterion will # Test completion
    #       terminated the given thread                                                                       # Test completion
    #       The criteria have priority among each other, which is the above specified order                   # Test completion
    #   If not done then replace the old bait file with the new bait                                          # Test completion
    # Inputs:                                                                                                 # Test completion
    my ($log_ref,       # The reference for the log file handle                                               # Test completion
        $ref_ref,       # The reference to the reference fasta file                                           # Test completion
        $old_bait_ref,  # The reference to the current bait file                                              # Test completion
        $type_ref,      # The type of analysis                                                                # Test completion
        $new_bait_ref,  # The new assembly file (future bait)                                                 # Test completion
        $min_len_ref,   # The minimum length for contigs to be kept                                           # Test completion
        $tab,           # Reference to the variable possibly holding some tabs                                # Test completion
        $asmbls_ref,    # Reference to array holding the assembly files that need to be modified if min_len   # Test completion
                        ) = @_;                                                                               # Test completion
                                                                                                              # Test completion
    # If $tab was not defined then use a reference to an empty string                                         # Test completion
    $tab = \"" unless $tab;#"                                                                                 # Test completion
                                                                                                              # Test completion
    # 1) Done IF there is NO bait file OR it is EMPTY                                                         # Test completion: No assembly
    unless (-e $$new_bait_ref) {                                                                              # Test completion: No assembly
        $print = strftime('%H:%M:%S',localtime) . "\t$$tab\tThere is no assembly\n";                          # Test completion: No assembly
        print {$$log_ref} $print; print $print;                                                               # Test completion: No assembly
        return "No assembly";                                                                                 # Test completion: No assembly
    }                                                                                                         # Test completion: No assembly
    if (-z $$new_bait_ref) {                                                                                  # Test completion: No assembly
        $print = strftime('%H:%M:%S',localtime) . "\t$$tab\tThere is no assembly\n";                          # Test completion: No assembly
        print {$$log_ref} $print; print $print;                                                               # Test completion: No assembly
        return "No assembly";                                                                                 # Test completion: No assembly
    }                                                                                                         # Test completion
                                                                                                              # Test completion
    # 2) Done if the old and new bait files are the same                                                      # Test completion: Same bait
    # Read in the bait files into hash variables                                                              # Test completion: Same bait
    my (%new_fasta, %old_fasta);                                                                              # Test completion: Same bait
    &fasta2hash($old_bait_ref, \%old_fasta, []);                                                              # Test completion: Same bait
    &fasta2hash($new_bait_ref, \%new_fasta, []);                                                              # Test completion: Same bait
                                                                                                              # Test completion: Same bait
    # Store the complete assembly for the statistics                                                          # Test completion: Same bait
    my %assembly = %new_fasta;                                                                                # Test completion: Same bait
    # A hash to store the fasta data for the new bait file                                                    # Test completion: Same bait
    my %bait_fasta;                                                                                           # Test completion: Same bait
                                                                                                              # Test completion:
    my %screened;                                                                                             # Test completion: Min length
    # Remove short contigs from the new one, old one was checked during the previous round                    # Test completion: Min length
    if ($$min_len_ref) {                                                                                      # Test completion: Min length
        for my $contig (keys %new_fasta) {                                                                    # Test completion: Min length
            if ($$min_len_ref < length $new_fasta{$contig}) {                                                 # Test completion: Min length
                $screened{$contig} = $new_fasta{$contig};                                                     # Test completion: Min length
            }                                                                                                 # Test completion: Min length
        }                                                                                                     # Test completion: Min length
        %new_fasta = %screened;                                                                               # Test completion: Min length
        for (@$asmbls_ref) {                                                                                  # Test completion: Min length
            open my $out_handle, '>', $_;                                                                     # Test completion: Min length
            for (keys %screened) {                                                                            # Test completion: Min length
                print {$out_handle} &to_fasta($_, $screened{$_});                                             # Test completion: Min length
            }                                                                                                 # Test completion: Min length
            close $out_handle;                                                                                # Test completion: Min length
        }                                                                                                     # Test completion: Same bait
    }                                                                                                         # Test completion: Same bait
    # The short contigs should be removed from both hashes                                                    # Test completion: Same bait
    %bait_fasta = %new_fasta;                                                                                 # Test completion: Same bait
                                                                                                              # Test completion: Same bait
    # Compare every contig of the bait and the new assembly with each other                                   # Test completion: Same bait
    #   If the sequence of a new contig is the same as that of a bait contig (or its reverse complement)      # Test completion: Same bait
    #       Then remove that contig from the hash                                                             # Test completion: Same bait
    #       If the hash is empty by the end then it means there is nothing new in the assembly                # Test completion: Same bait
    for my $a (keys %new_fasta) {                                                                             # Test completion: Same bait
        for my $b (keys %old_fasta) {                                                                         # Test completion: Same bait
            # if a contig has already been removed then skip to the next one                                  # Test completion: Same bait
            next unless $new_fasta{$a};                                                                       # Test completion: Same bait
            my $reverse = reverse $old_fasta{$b};                                                             # Test completion: Same bait
            $reverse =~ tr/ATCG/TAGC/;                                                                        # Test completion: Same bait
            if ($new_fasta{$a} eq $old_fasta{$b}) {                                                           # Test completion: Same bait
                $new_fasta{$a} = undef;                                                                       # Test completion: Same bait
            } elsif ($new_fasta{$a} eq $reverse) {                                                            # Test completion: Same bait
                $new_fasta{$a} = undef;                                                                       # Test completion: Same bait
            }                                                                                                 # Test completion: Same bait
        }                                                                                                     # Test completion: Same bait
    }                                                                                                         # Test completion: Same bait
                                                                                                              # Test completion: Same bait
                                                                                                              # Test completion
    # Print some information about the assembly                                                               # Test completion
    my $statistics = &assembly_statistics(\%assembly, $tab);                                                  # Test completion
    $print = $statistics;                                                                                     # Test completion
    print {$$log_ref} $print; print $print;                                                                   # Test completion
                                                                                                              # Test completion
    if ($$min_len_ref) {                                                                                      # Test completion
        my $statistics = &assembly_statistics(\%screened, $tab);                                              # Test completion
        $print = strftime('%H:%M:%S',localtime);                                                              # Test completion
        $print = "\t$$tab\tAssembly statistics after the removal of short contigs\n";                         # Test completion
        $print .= $statistics;                                                                                # Test completion
        print {$$log_ref} $print; print $print;                                                               # Test completion
                                                                                                              # Test completion
                                                                                                              # Test completion
        # unless there is a value that is defined, there was nothing new in the assembly so stop the thread   # Test completion: Same bait
        unless (grep{$_} values %screened) {                                                                  # Test completion: Same bait
            $print = strftime('%H:%M:%S',localtime);                                                          # Test completion: Same bait
            $print = "\t$$tab\tAssembly contains no contigs that are of sufficient length\n";                 # Test completion: Same bait
            print {$$log_ref} $print; print $print;                                                           # Test completion: Same bait
            return "Same bait";                                                                               # Test completion: Same bait
        }                                                                                                     # Test completion: Same bait
    }                                                                                                         # Test completion: Same bait
                                                                                                              # Test completion: Same bait
    # unless there is a value that is defined, there was nothing new in the assembly so stop the thread       # Test completion: Same bait
    unless (grep{$_} values %new_fasta) {                                                                     # Test completion: Same bait
        $print = strftime('%H:%M:%S',localtime) . "\t$$tab\tAssembly is the same as the bait\n";              # Test completion: Same bait
        print {$$log_ref} $print; print $print;                                                               # Test completion: Same bait
        return "Same bait";                                                                                   # Test completion: Same bait
    }                                                                                                         # Test completion
                                                                                                              # Test completion
    # 3) Done if homologous region is completely matched using exonerate or there is no improvement           # Test completion: Exonerate
    if ($$type_ref && $$type_ref =~ /exonerate/) {                                                            # Test completion: Exonerate
        $print = strftime('%H:%M:%S',localtime) . "\t$$tab\tCheck completion using exonerate\n";              # Test completion: Exonerate
        print {$$log_ref} $print; print $print;                                                               # Test completion: Exonerate
                                                                                                              # Test completion: Exonerate
        # Run the exonerate subroutine and print the result if it returns true                                # Test completion: Exonerate
        my $test = &exonerate($ref_ref, $new_bait_ref, $tab, $log_ref);                                       # Test completion: Exonerate
        if ($test) {                                                                                          # Test completion: Exonerate
            my $time = strftime('%H:%M:%S',localtime);                                                        # Test completion: Exonerate
            if ($test eq "No improvement") {                                                                  # Test completion: Exonerate
                $print = "$time\t$$tab\tThere was no improvement, this thread is terminated\n";               # Test completion: Exonerate
                print {$$log_ref} $print; print $print;                                                       # Test completion: Exonerate
            } elsif ($test eq "Multiple contigs") {                                                           # Test completion: Exonerate
                $print = "$time\t$$tab\tThe reference matched to multiple contigs of the assembly\n";         # Test completion: Exonerate
                print {$$log_ref} $print; print $print;                                                       # Test completion: Exonerate
            } else {                                                                                          # Test completion: Exonerate
                $print = "$time\t$$tab\tThe reference matched to $test in the assembly\n";                    # Test completion: Exonerate
                print {$$log_ref} $print; print $print;                                                       # Test completion: Exonerate
                $test = "Found";                                                                              # Test completion: Exonerate
            }                                                                                                 # Test completion: Exonerate
            return $test;                                                                                     # Test completion: Exonerate
        }                                                                                                     # Test completion: Exonerate
    }                                                                                                         # Test completion
                                                                                                              # Test completion
    # 4) Done if the total assembly is at least a given length                                                # Test completion: Length (total)
    if ($$type_ref && $$type_ref =~ /total=(\d+)/i) {                                                         # Test completion: Length (total)
        my $min = $1;                                                                                         # Test completion: Length (total)
        my $length_spec = "total length is at least $min bp";                                                 # Test completion: Length (total)
        $statistics =~ /Total size is:\s+(\d+)/;                                                              # Test completion: Length (total)
        my $value = $1;                                                                                       # Test completion: Length (total)
        if ($value >= $min) {                                                                                 # Test completion: Length (total)
            $print = "\t$$tab\tThe assembly satisfied the length criterion: $length_spec\n";                  # Test completion: Length (total)
            print {$$log_ref} $print; print $print;                                                           # Test completion: Length (total)
            return "Reached: $length_spec";                                                                   # Test completion: Length (total)
        }                                                                                                     # Test completion: Length (total)
    }                                                                                                         # Test completion: Length
                                                                                                              # Test completion
    # 5) Done if the longest contig is at least a given length                                                # Test completion: Length (longest)
    if ($$type_ref && $$type_ref =~ /longest=(\d+)/i) {                                                       # Test completion: Length (longest)
        my $min = $1;                                                                                         # Test completion: Length (longest)
        my $length_spec = "longest contig is at least $min bp";                                               # Test completion: Length (longest)
        $statistics =~ /Longest contig is:\s+(\d+)/;                                                          # Test completion: Length (longest)
        my $value = $1;                                                                                       # Test completion: Length (longest)
        if ($value >= $min) {                                                                                 # Test completion: Length (longest)
            $print = "\t$$tab\tThe assembly satisfied the length criterion: $length_spec\n";                  # Test completion: Length (longest)
            print {$$log_ref} $print; print $print;                                                           # Test completion: Length (longest)
            return "Reached: $length_spec";                                                                   # Test completion: Length (longest)
        }                                                                                                     # Test completion: Length (longest)
    }                                                                                                         # Test completion: Length
                                                                                                              # Test completion
    # 6) Done if the N50 value of the assembly is at least a given length                                     # Test completion: Length (N50)
    if ($$type_ref && $$type_ref =~ /N50=(\d+)/i) {                                                           # Test completion: Length (N50)
        my $min = $1;                                                                                         # Test completion: Length (N50)
        my $length_spec = "N50 is at least $min bp";                                                          # Test completion: Length (N50)
        $statistics =~ /N50 is:\s+(\d+)/;                                                                     # Test completion: Length (N50)
        my $value = $1;                                                                                       # Test completion: Length (N50)
        if ($value >= $min) {                                                                                 # Test completion: Length (N50)
            $print = "\t$$tab\tThe assembly satisfied the length criterion: $length_spec\n";                  # Test completion: Length (N50)
            print {$$log_ref} $print; print $print;                                                           # Test completion: Length (N50)
            return "Reached: $length_spec";                                                                   # Test completion: Length (N50)
        }                                                                                                     # Test completion: Length
    }                                                                                                         # Test completion
                                                                                                              # Test completion
    # No criteria was met so copy the assembly file to be the new bait file                                   # Test completion: continue
    open my $bait_handle, '>', $$old_bait_ref;                                                                # Test completion: continue
    for (keys %bait_fasta) {                                                                                  # Test completion: continue
        print {$bait_handle} &to_fasta($_, $bait_fasta{$_});                                                  # Test completion: continue
    }                                                                                                         # Test completion: continue
    close $bait_handle;                                                                                       # Test completion: continue
    return;                                                                                                   # Test completion: continue
}                                                                                                             # Test completion
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
sub assembly_statistics{                                                                                      # Test completion: statistics
    # This subroutine returns a string containing some statistics on the input assembly                       # Test completion: statistics
    #   Inputs:                                                                                               # Test completion: statistics
    my ($assembly_ref,	# The reference to the assembly file                                                  # Test completion: statistics
        $tab            # Reference to the variable possibly holding some tabs                                # Test completion: statistics
                        ) = @_;                                                                               # Test completion: statistics
                                                                                                              # Test completion: statistics
    # If $tab was not defined then use a reference to an empty string                                         # Test completion: statistics
    $tab = \"" unless $tab;#"                                                                                 # Test completion: statistics
                                                                                                              # Test completion: statistics
    # A variable to store the return string                                                                   # Test completion: statistics
    my $print;                                                                                                # Test completion: statistics
                                                                                                              # Test completion: statistics
    # Get the number of contigs in the assembly                                                               # Test completion: statistics
    my $contig_count = keys %$assembly_ref;                                                                   # Test completion: statistics
                                                                                                              # Test completion: statistics
    # sort the contigs in descending order based on contig size                                               # Test completion: statistics
    my @descending = sort{$b<=>$a} map{length $_} values %$assembly_ref;                                      # Test completion: statistics
                                                                                                              # Test completion: statistics
    # Get the total length of the assembly                                                                    # Test completion: statistics
    my $total;                                                                                                # Test completion: statistics
    for (@descending) {                                                                                       # Test completion: statistics
        $total += $_;                                                                                         # Test completion: statistics
    }                                                                                                         # Test completion: statistics
                                                                                                              # Test completion: statistics
    # Get the N50 of the assembly                                                                             # Test completion: statistics
    my $n50;                                                                                                  # Test completion: statistics
    my $accumulator;                                                                                          # Test completion: statistics
    my $half = $total / 2;                                                                                    # Test completion: statistics
    for (@descending) {                                                                                       # Test completion: statistics
        $accumulator += $_;                                                                                   # Test completion: statistics
        if ($accumulator >= $half) {                                                                          # Test completion: statistics
            $n50 = $_;                                                                                        # Test completion: statistics
            last;                                                                                             # Test completion: statistics
        }                                                                                                     # Test completion: statistics
    }                                                                                                         # Test completion: statistics
                                                                                                              # Test completion: statistics
    # Get the mean contig length                                                                              # Test completion: statistics
    my $mean = $total / $contig_count;                                                                        # Test completion: statistics
                                                                                                              # Test completion: statistics
    # Get the shortest and longest contig lengths                                                             # Test completion: statistics
    my ($min, $max) = (@descending)[-1,0];                                                                    # Test completion: statistics
                                                                                                              # Test completion: statistics
    # Save the data to the return string                                                                      # Test completion: statistics
    $print .= strftime('%H:%M:%S',localtime) . "\t$$tab\tAssembly statistics (in bp):\n";                     # Test completion: statistics
    $print .= strftime('%H:%M:%S',localtime) . "\t$$tab\t\tNumber of contigs:\t$contig_count\n";              # Test completion: statistics
    $print .= strftime('%H:%M:%S',localtime) . "\t$$tab\t\tTotal size is:\t\t$total\n";                       # Test completion: statistics
    $print .= strftime('%H:%M:%S',localtime) . "\t$$tab\t\tLongest contig is:\t$max\n";                       # Test completion: statistics
    $print .= strftime('%H:%M:%S',localtime) . "\t$$tab\t\tN50 is:\t\t\t$n50\n";                              # Test completion: statistics
    $print .= strftime('%H:%M:%S',localtime) . "\t$$tab\t\tMean contig size is:\t$mean\n";                    # Test completion: statistics
    $print .= strftime('%H:%M:%S',localtime) . "\t$$tab\t\tShortest contig is:\t$min\n";                      # Test completion: statistics
    return $print;                                                                                            # Test completion: statistics
}                                                                                                             # Test completion: statistics
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
sub exonerate{                                                                                                # Test completion: Exonerate
    # This subroutine uses exonerate to align the reference and the assembly to see if the target is fully    # Test completion: Exonerate
    #   recovered. Can have 4 different return values.                                                        # Test completion: Exonerate
    #       False if there was improvement in the matching                                                    # Test completion: Exonerate
    #       Name of the contig that matched the reference completely (the matched region is saved to file)    # Test completion: Exonerate
    #       "Multiple contigs" the reference matched to multiple contigs, then it could possibly recovered    # Test completion: Exonerate
    #       "No improvement" if the same positions are still unmatched                                        # Test completion: Exonerate
    # The output of exonerate run is saved to "$$q_ref\.exonerate"                                            # Test completion: Exonerate
    #   Inputs:                                                                                               # Test completion: Exonerate
    my ($q_ref,     # The reference to the reference fasta file                                               # Test completion: Exonerate
        $t_ref,     # The reference to the assembly fasta file                                                # Test completion: Exonerate
        $tab,       # Reference to the variable possibly holding some tabs                                    # Test completion: Exonerate
        $log_ref,   # The reference for the log file handle                                                   # Test completion: Exonerate
                    ) = @_;                                                                                   # Test completion: Exonerate
                                                                                                              # Test completion: Exonerate
    # If $tab was not defined then use a reference to an empty string                                         # Test completion: Exonerate
    $tab = \"" unless $tab;#"                                                                                 # Test completion: Exonerate
                                                                                                              # Test completion: Exonerate
    # A variable to store the return string                                                                   # Test completion: Exonerate
    my $print;                                                                                                # Test completion: Exonerate
                                                                                                              # Test completion: Exonerate
    # The string to store witch positions are matched (+) and unmatched (-) in the reference                  # Test completion: Exonerate
    my $matched;                                                                                              # Test completion: Exonerate
                                                                                                              # Test completion: Exonerate
    # The matched string from the previous round                                                              # Test completion: Exonerate
    my $matched_old = "";                                                                                     # Test completion: Exonerate
    if (-e "$$q_ref\.match") {                                                                                # Test completion: Exonerate
        open my $previous, '<', "$$q_ref\.match";                                                             # Test completion: Exonerate
        $matched_old = <$previous>;                                                                           # Test completion: Exonerate
        close $previous;                                                                                      # Test completion: Exonerate
    }                                                                                                         # Test completion: Exonerate
                                                                                                              # Test completion: Exonerate
    # Command to run the exonerate analysis                                                                   # Test completion: Exonerate
    my $cmd = "$exonerate_cmd $$q_ref $$t_ref --percent 5  --model affine:overlap -E --ryo ";                 # Test completion: Exonerate
    $cmd .= '"(%qi %ql %qab %qae %qal %qS) (%ti %tl %tab %tae %tal %tS) (%V) (%et %es %em %pi)\n" ';          # Test completion: Exonerate
    $cmd .= "--showalignment no --showvulgar no";                                                             # Test completion: Exonerate
                                                                                                              # Test completion: Exonerate
    # Run the command                                                                                         # Test completion: Exonerate
    my $exonerate = `$cmd 2>exonerate.log`;                                                                   # Test completion: Exonerate
                                                                                                              # Test completion: Exonerate
    # Save the exonerate output to file                                                                       # Test completion: Exonerate
    open my $exonerate_save, '>', "$$q_ref\.exonerate";                                                       # Test completion: Exonerate
    print {$exonerate_save} $exonerate;                                                                       # Test completion: Exonerate
    close $exonerate_save;                                                                                    # Test completion: Exonerate
                                                                                                              # Test completion: Exonerate
    # Parse the exonerate result                                                                              # Test completion: Exonerate
    for (split /\n/, $exonerate) {                                                                            # Test completion: Exonerate
        if (/\((\S+)\ (\d+)\ (\d+)\ (\d+)\ (\d+)\ ([-+])\)\ \((\S+)\ (\d+)\ (\d+)\ (\d+)\ (\d+)\ ([-+])\)/x){ # Test completion: Exonerate
           # ^ query (ID length start end aligned strand)   ^ target (ID length start end aligned strand)     # Test completion: Exonerate
            my %query = (     id      => $1,                                                                  # Test completion: Exonerate
                              len     => $2,                                                                  # Test completion: Exonerate
                              start   => $3,                                                                  # Test completion: Exonerate
                              end     => $4,                                                                  # Test completion: Exonerate
                              aligned => $5,                                                                  # Test completion: Exonerate
                              strand  => $6,                                                                  # Test completion: Exonerate
                                              );                                                              # Test completion: Exonerate
            my %target = (    id      => $7,                                                                  # Test completion: Exonerate
                              len     => $8,                                                                  # Test completion: Exonerate
                              start   => $9,                                                                  # Test completion: Exonerate
                              end     => $10,                                                                 # Test completion: Exonerate
                              aligned => $11,                                                                 # Test completion: Exonerate
                              strand  => $12,                                                                 # Test completion: Exonerate
                                              );                                                              # Test completion: Exonerate
                                                                                                              # Test completion: Exonerate
            # Make the query positive stranded                                                                # Test completion: Exonerate
            #   If needed we have to reverse the target information and strand                                # Test completion: Exonerate
            if ($query{"strand"} eq "-") {                                                                    # Test completion: Exonerate
                ($query{"start"}, $query{"end"}) = ($query{"end"}, $query{"start"});                          # Test completion: Exonerate
                $query{"strand"} = "+";                                                                       # Test completion: Exonerate
                ($target{"start"}, $target{"end"}) = ($target{"end"}, $target{"start"});                      # Test completion: Exonerate
                if ($target{"strand"} eq "-") {                                                               # Test completion: Exonerate
                    $target{"strand"}  = "+";                                                                 # Test completion: Exonerate
                } else {                                                                                      # Test completion: Exonerate
                    $target{"strand"}  = "-";                                                                 # Test completion: Exonerate
                }                                                                                             # Test completion: Exonerate
            }                                                                                                 # Test completion: Exonerate
                                                                                                              # Test completion: Exonerate
            # Initialize the matched string with no bases matched                                             # Test completion: Exonerate
            $matched = "-" x $query{len};                                                                     # Test completion: Exonerate
                                                                                                              # Test completion: Exonerate
            # Did the current hit completely cover the reference                                              # Test completion: Exonerate
            if ( ($query{start} == 0 && $query{len} == $query{end}) ||                                        # Test completion: Exonerate
                 ($query{end} == 0 && $query{len} == $query{start})    ) {                                    # Test completion: Exonerate
                # The query is completely matched to one of the contigs                                       # Test completion: Exonerate
                                                                                                              # Test completion: Exonerate
                # Create the result file                                                                      # Test completion: Exonerate
                #   Read in the assembly file                                                                 # Test completion: Exonerate
                my %assembly;                                                                                 # Test completion: Exonerate
                &fasta2hash($t_ref, \%assembly, []);                                                          # Test completion: Exonerate
                                                                                                              # Test completion: Exonerate
                # Get the sequence of the contig that matched                                                 # Test completion: Exonerate
                #   Reverse complement it if the negative strand was matched                                  # Test completion: Exonerate
                my $seq;                                                                                      # Test completion: Exonerate
                for (keys %assembly) {                                                                        # Test completion: Exonerate
                    if ($_ =~ /^$target{id}\b/) {                                                             # Test completion: Exonerate
                        if ($target{"strand"} eq "+") {                                                       # Test completion: Exonerate
                            $seq = substr($assembly{$_}, $target{"start"}, $target{"aligned"});               # Test completion: Exonerate
                        } else {                                                                              # Test completion: Exonerate
                            $seq = reverse substr($assembly{$_}, $target{"end"}, $target{"aligned"});         # Test completion: Exonerate
                            $seq =~ tr/ACTG/TGAC/;                                                            # Test completion: Exonerate
                        }                                                                                     # Test completion: Exonerate
                        last;                                                                                 # Test completion: Exonerate
                    }                                                                                         # Test completion: Exonerate
                }                                                                                             # Test completion: Exonerate
                                                                                                              # Test completion: Exonerate
                # Save the matched region to file                                                             # Test completion: Exonerate
                open my $result, '>', "result.fas";                                                           # Test completion: Exonerate
                print {$result} &to_fasta($target{"id"}, $seq);                                               # Test completion: Exonerate
                close $result;                                                                                # Test completion: Exonerate
                                                                                                              # Test completion: Exonerate
                # Return the name of the contig                                                               # Test completion: Exonerate
                return $target{"id"};                                                                         # Test completion: Exonerate
            } else {                                                                                          # Test completion: Exonerate
                                                                                                              # Test completion: Exonerate
                # Get the start and end of the matched region                                                 # Test completion: Exonerate
                my ($a, $b);                                                                                  # Test completion: Exonerate
                if ($query{strand} eq "+") {                                                                  # Test completion: Exonerate
                    $a = $query{start};                                                                       # Test completion: Exonerate
                    $b = $query{end};                                                                         # Test completion: Exonerate
                } else {                                                                                      # Test completion: Exonerate
                    $a = $query{end};                                                                         # Test completion: Exonerate
                    $b = $query{start};                                                                       # Test completion: Exonerate
                }                                                                                             # Test completion: Exonerate
                                                                                                              # Test completion: Exonerate
                # Store the new matched string in a separate variable for a few lines                         # Test completion: Exonerate
                my $matched_new;                                                                              # Test completion: Exonerate
                                                                                                              # Test completion: Exonerate
                # Copy the matching information for the bases that weren't matched by this contig             # Test completion: Exonerate
                if ($a > 0) {                                                                                 # Test completion: Exonerate
                    $matched_new = substr($matched, 0, $a);                                                   # Test completion: Exonerate
                }                                                                                             # Test completion: Exonerate
                                                                                                              # Test completion: Exonerate
                # Mark the bases as matched (+) for the bases that have been matched by this contig           # Test completion: Exonerate
                $matched_new .= "+" x ($b - $a);                                                              # Test completion: Exonerate
                                                                                                              # Test completion: Exonerate
                # Copy the matching information for the bases that weren't matched by this contig             # Test completion: Exonerate
                if ($b < $query{len}) {                                                                       # Test completion: Exonerate
                    $matched_new .= substr($matched, $b);                                                     # Test completion: Exonerate
                }                                                                                             # Test completion: Exonerate
                $matched = $matched_new;                                                                      # Test completion: Exonerate
            }                                                                                                 # Test completion: Exonerate
        }                                                                                                     # Test completion: Exonerate
    }                                                                                                         # Test completion: Exonerate
                                                                                                              # Test completion: Exonerate
    # If it got this far there was no single good hit, but still the reference is completely matched          # Test completion: Exonerate
    #   That means that the reference matched to multiple contigs                                             # Test completion: Exonerate
    if ($matched !~ /-/) {                                                                                    # Test completion: Exonerate
        # All the positions were matched                                                                      # Test completion: Exonerate
        return "Multiple contigs";                                                                            # Test completion: Exonerate
    }                                                                                                         # Test completion: Exonerate
                                                                                                              # Test completion: Exonerate
    # Check if there is some change in the matched string from last round                                     # Test completion: Exonerate
    if ($matched eq $matched_old) {                                                                           # Test completion: Exonerate
        # If no then there is nothing new that helps us                                                       # Test completion: Exonerate
        return "No improvement";                                                                              # Test completion: Exonerate
    } else {                                                                                                  # Test completion: Exonerate
        # There is still some hope, because there was some change                                             # Test completion: Exonerate
        open my $save, '>', "$$q_ref\.match";                                                                 # Test completion: Exonerate
        print {$save} "$matched";                                                                             # Test completion: Exonerate
        close $save;                                                                                          # Test completion: Exonerate
        # Print some information about the match                                                              # Test completion: Exonerate
	my $ref_len = length $matched;                                                                        # Test completion: Exonerate
	my $ref_m = length join "", grep {$_ eq "+"} split //, $matched;                                      # Test completion: Exonerate
        $print .= strftime('%H:%M:%S',localtime) . "\t$$tab\t\t$ref_m bp are matched out of $ref_len bp\n";   # Test completion: Exonerate
        print {$$log_ref} $print; print $print;                                                               # Test completion: Exonerate        
    }                                                                                                         # Test completion: Exonerate
                                                                                                              # Test completion: Exonerate
    return undef;                                                                                             # Test completion: Exonerate
}                                                                                                             # Test completion: Exonerate
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
sub fasta2hash{                                                                                               # fasta2hash
    # This subroutine reads in a fasta file (input) and saves it into a hash (output), saves ids to an array  # fasta2hash
    #   keys: ids    values: sequences                                                                        # fasta2hash
    #   Inputs:                                                                                               # fasta2hash
    my ($file_ref,  # Reference to the FASTA file                                                             # fasta2hash
        $hash_ref,  # Reference to the hash variable that will hold the entries of the Fasta file             # fasta2hash
                    #   keys: id line       values: the corresponding sequence                                # fasta2hash
        $order_ref  # Reference to an array that stores the order of the entries (elements are the id lines)  # fasta2hash
                    ) = @_;                                                                                   # fasta2hash
                                                                                                              # fasta2hash
    # A scalar variable to store the id of the current entry                                                  # fasta2hash
    my $current;                                                                                              # fasta2hash
                                                                                                              # fasta2hash
    # Open the FASTA file for reading, then read the file in line by line                                     # fasta2hash
    open my $fasta, '<', $$file_ref;                                                                          # fasta2hash
    for (<$fasta>) {                                                                                          # fasta2hash
        s/\R//g;                        # Remove line endings                                                 # fasta2hash
        next if /^\s*$/;                # Skip empty lines                                                    # fasta2hash
        if (/^>(.*)/) {                 # Get the ID                                                          # fasta2hash
            $current = $1;              # save the ID                                                         # fasta2hash
            push @$order_ref, $current;	# save the ID in the order array                                      # fasta2hash
        } else {                                                                                              # fasta2hash
            $$hash_ref{$current} .= $_;	# append line to the current sequence                                 # fasta2hash
        }                                                                                                     # fasta2hash
    }                                                                                                         # fasta2hash
    close $fasta;                                                                                             # fasta2hash
}                                                                                                             # fasta2hash
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#======================================END===================================================================================================
# The MIT License (MIT)
#
# Copyright (c) 2015 b-brankovics
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
