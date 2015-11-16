# GRAbB
GRAbB (Genome Region Assembly by Baiting) is program designed to assemble selected regions of the genome or transcriptome using reference sequences and NGS data.

-------------------------------------------
# Table of contents


1. [Usage](#usage)
    1. [Installation](#installation)
    2. [Documentation](#documentation)
    3. [Examples](#examples)
2. [Prerequisites](#prerequisites)
    1. [mirabait](#mirabait)
    2. [Edena](#edena)
    3. [Velvet](#velvet)
    4. [Seqtk](#seqtk)
    5. [exonerate](#exonerate)
    6. [PRINSEQ lite](#prinseq-lite)
3. [Helper programs](#helper-programs)
    1. [fastq2fasta](#fastq2fasta)
    2. [get_overlaps](#get_overlaps)
    3. [interleaved2pairs](#interleaved2pairs)
    4. [merge_contigs](#merge_contigs)
    5. [rename_fastq](#rename_fastq)
    6. [single2pairs](#single2pairs)
    7. [uniform_length](#uniform_length)
    8. [fasta_shift](#fasta_shift)
    9. [pairwise_alignment_score](#pairwise_alignment_score)
    10. [reverse_complement](#reverse_complement)
4. [Algorithm overview](#algorithm-overview)
    1. [Main loop](#main-loop)
    2. [Creating the bait](#creating-the-bait)
    3. [Baiting](#baiting)
    4. [Collecting reads](#collecting-reads)
    5. [Assembly](#assembly)
    6. [Testing completion](#testing-completion)
    7. [Modes](#modes)
5. [Arguments](#arguments)
    1. [ref](#ref)
    2. [bait](#bait)
    3. [reads](#reads)
    4. [folder](#folder)
    5. [prefix](#prefix)
    6. [single](#single)
    7. [min_length](#min_length)
    8. [type](#type)
    9. [arg1](#arg1)
    10. [arg2](#arg2)
    11. [assembler](#assembler)
    12. [clean](#clean)
6. [Using custom assembler program](#using-custom-assembler-program)
    1. [Adding to the source code of GRAbB](#adding-to-the-source-code-of-grabb)
    2. [Using external_skeleton](#using-external_skeleton)
7. [Citation](#citation)
8. [Contact](#contact)

----------------------------
#### Usage


#### Installation

1. Install prerequisites (If this step is skipped, then
   **configure_GRAbB.pl** tries to use prerequisites included in the package)
    - Minimal set:
        + Baiting program: **mirabait** _(recommended)_ OR
		kmer\_bait\.pl (no installation needed)
	    + Read collecting program: **seqtk** _(recommended)_ OR
          create\_readpool\.pl (no installation needed)
		+ Assembler: Edena OR Vevlet OR external\_scaffold\.pl (needs
          to be modified and requires a working installation of an assembler)

    - Assemblers:
		+ Edena: default assembler for GRAbB.pl
		+ Vevlet
		+ Other assembler: external\_scaffold\.pl has to be also edited

    - Exonerate:
		+ Required for running GRAbB.pl in exonerate mode

2. Configure GRAbB.pl

        ./configure_GRAbB.pl

    Recommended to add prerequisites to the Path Or add the absolute
    path of the executables to the source code of GRAbB.pl before
	running configure_GRAbB.pl

    Configured GRAbB.pl can be found in bin directory

    **Bug:** On some systems the exonerate binary included in the
      package runs extremely slow. Configure can get stuck at 'Testing
      exonerate'
	  block, then issue **Ctrl + c**.
	  On these systems exonerate has to be installed or build from
      source code (See [Prerequisites::exonerate](#exonerate), before
      rerunning the configuration script.

2. To test installation run the following (assembler has to be
   adjusted unless GRAbB.pl is configured without **Edena**)

        bin/GRAbB.pl --ref for_testing/assembly.fas --reads for_testing/read* --folder test --prefix test

OR

Use Docker (See __[Docker.md](Docker.md)__ for more detailed instructions)

* Either download the docker repository via docker interface `docker pull brankovics/grabb`
* Or create a local docker image:

        git clone https://github.com/b-brankovics/grabb
        cd grabb/docker
        docker build -t localhost:5000/$USER/grabb .

#### Documentation

Run GRAbB.pl without any arguments and it prints the Usage information

The documentation is this file and the files mentioned at the [examples](#examples).

#### Examples

See the wiki or the files __[Docker.md](Docker.md)__, __[Examples.md](Examples.md)__ and __[Tutorial.md](Tutorial.md)__.

--------------------------------
#### Prerequisites

#### mirabait

1. Download MIRA (4.0) assembler (http://sourceforge.net/projects/mira-assembler/)
2. Extract it. The executable files are in the _bin_ folder.
3. Copy/move or symlink '**mirabait**' into somewhere in the path or add the folder to the path (This program uses only mirabait)

    __Warning__: the name of the executable file has to be mirabait!

    __Bug__: mira 4.0.2 does not work properly, but mira 4.0 does

OR you may also use [kmer\_bait\.pl](perl_programs/kmer_bait.pl),
which is less efficient, but uses only perl and standard Unix commands

#### Edena

1. Download [EDENA](http://www.genomic.ch/edena.php) and extract it or use the copy in the 3rd_party_programs
2. Change to the directory
3. Type `make` on the command line (g++ needs to be installed, on ubuntu type `sudo apt-get install g++`)
4. Copy/move or symlink '**edena**' into somewhere in the path or add the folder to the path (The files in the bin folder)

#### Velvet

1. Download [Velvet](https://www.ebi.ac.uk/~zerbino/velvet/) and extract it or use the copy in the 3rd_party_programs
2. Change to the directory. First zlib needs to be installed
3. Change to 'third-party/zlib-1.2.3/'
4. Type `make` on the command line
5. Type `sudo make install` on the command line
6. Go back to the parent directory (`cd ../..`)
7. Type `make` on the command line
8. Copy/move or symlink '__velveth__' and '__velvetg__' into somewhere in the path or add the folder to the path

#### Seqtk

1. Download [Seqtk](https://github.com/lh3/seqtk.git) from github and uncompress it or `git clone https://github.com/lh3/seqtk.git`
2. Change to the directory
3. Type `make` on the command line (zlib needs to be installed, see 1.3) for instructions)
4. Copy/move or symlink '__seqtk__' into somewhere in the path or add the folder to the path (The files in the bin folder)

OR you may also use [create\_readpool\.pl](perl_programs/create_readpool.pl),
which is less efficient, but uses only perl and standard Unix commands

#### exonerate

For Ubuntu run `sudo apt-get install exonerate`

Else:

1. Download [exonerate](https://www.ebi.ac.uk/sites/ebi.ac.uk/files/groups/flicek/exonerate/exonerate-2.2.0.tar.gz)
   from the EBI website and uncompress it or use the version included
   in the GRAbB package (3rd\_party\_programs)
2. Change to the directory
3. Type the following commands (The following packages have to be
   installed on the system before running `./configure`: **gcc**, **make** and **glib2**)

	    ./configure
	    make
	    make check
	    make install
4. The executable is found at `src/program/exonerate`. Copy/move or
   symlink '__exonerate__' into somewhere in the path or add the folder to
   the path (The files in the src/program folder)

#### PRINSEQ lite

1. Download [PRINSEQ lite](http://prinseq.sourceforge.net/) and extract it or use the copy in the 3rd_party_programs
2. Copy/move or symlink '**prinseq-lite.pl**' into somewhere in the path or add the folder to the path

-----------------------
#### Helper programs

#### fastq2fasta

This program creates a FASTA format read file for each FASTQ read file specified.

Usage:

    ./fastq2fasta <reads_1.fastq> <reads_2.fastq>

#### get_overlaps

This program reads the contigs from a fasta file and checks if they
are overlapping with each other by using a minimal overlap size that is
specified at invocation. Finally, prints all the overlaps found.

Usage:

    ./get_overlaps <contigs.fasta> <overlap>

#### interleaved2pairs

This program creates a forward and reverse read file from an interleaved file

Usage:

    ./interleaved2pairs <reads.fastq>

#### merge_contigs

This program reads the contigs from a fasta file and checks if they
are overlapping with each other by using a minimal overlap size that is
specified at invocation. Afterwards it loops through all the contigs
and merges contig pairs that only overlap with each other at the given
side. In the end it saves the contigs that were created into a file.

Usage:

    ./merge_contigs <contigs.fasta> <overlap> <output.fas>

#### rename_fastq

This program creates a FASTA format read file for each FASTQ read file

Usage:

    ./rename_fastq <reads_1.fastq> <reads_2.fastq>

#### single2pairs

This program creates paired-end read files from single-end files

1. Selects the reads that are at least as long as the specified length (\<_int_\>)
2. Gets the first \<_int_\> characters to be used as forward read
3. Gets the last \<_int_\> characters to be used as reverse read (reverse complement)

All the produced reads are $len long.

The output file will be created in the current working directory

Usage:

    ./single2pairs <int> <reads.fastq>

#### uniform_length

This program creates a read file with reads with uniform length

1. Selects the reads that are at least as long as the specified length (\<_int_\>)
2. Trims the reads to the specified length

All the produced reads are $len long

The output file will be created in the current working directory

Usage:

    ./uniform_length <int> <reads.fastq>

#### fasta_shift

This program takes a fasta file and shifts the sequence in it
according to a position value or a reference file the output is
printed to STDOUT

1. Using position value n is the position value

    Usage:

        ./fasta_shift -i <input.fas> -p <int> ><output.fas>

2. Using a reference file: the output will start with the first
occurrence of the sequence in the reference file

    Usage:

        ./fasta_shift -i <input.fas> -r <ref.fas> ><output.fas>

#### pairwise_alignment_score

This program compares two sequences and prints a few metrics

Output: <total length>   <# of identical bases>   <percentage of identity>  <# of identical bases, not counting '-'>  <percentage of identity w/o '-'>

Usage:

     ./pairwise_alignment_score <input1.fas> <input2.fas>

#### reverse_complement

This program prints the reverse complement of the input to STDOUT

The input may be files or STDIN

Usage:

     ./reverse_complement <input.fas> ><output.fas>

--------------------------------
#### Algorithm overview

GRAbB is written in Perl and it uses only modules that are part of the
core distribution of Perl. In addition to basic UNIX commands, the
following third-party programs are used by GRAbB: mirabait (from the
[MIRA](http://sourceforge.net/projects/mira-assembler/) package),
[Seqtk](https://github.com/lh3/seqtk),
[EDENA](http://www.genomic.ch/edena.php),
[Velvet](https://www.ebi.ac.uk/~zerbino/velvet/) and
[Exonerate](https://www.ebi.ac.uk/~guy/exonerate/).

The program is designed to be versatile and flexible with the
following functionalities:

+ Use pairing information
+ Use additional bait sequences
+ Assemble multiple regions separately in a single run
+ Use any of a range of completion criteria (may use different ones for each region)

These functionalities are detailed below at the appropriate steps of
the algorithm.

#### Main loop

The main loop of GRAbB can be summarized as follows:

1. creating the bait file
2. finding reads by baiting
3. collecting reads
4. _de novo_ assembly of selected reads
5. testing completion

If [multi-mode](#modes) is selected, then the general
baiting step is followed by specific baiting, _de novo_ assembling and
completion testing for each of the threads. The threads are generated
by splitting the reference file into single-entry FASTA files. These
newly created reference files are used as bait for the initial
specific baiting steps for the given thread. At the end of each cycle
of the main loop, the program checks whether there is any thread that
is not completed yet, then it continues or stops accordingly.

#### Creating the bait

At the invocation of GRAbB it is possible to specify a [length
filter](#min_length) that excludes all contigs from
the assembly that are shorter than the specified length from being
used for generating the bait file for the next iteration.

+ Initial bait file

    At the start of the run an initial (general) bait file is created
    by concatenating the [reference](#ref) and
    [bait](#bait) files. This bait files is used for
    the first general baiting step.

+ General bait file

    At the start of each iteration a general bait file is created by
    concatenating the latest assembly file(s)

+ Specific bait file

    The specific bait file is only created in
    [multi-mode](#modes). At the start of the run the
    reference file is split into single-entry FASTA files and these
    are used as initial specific bait files. In latter iterations the
    latest assembly for the give thread is used as the new bait file.

#### Baiting

Reads belonging to the specified sequence are identified by using
exact k-mer (31 bp) matching that is implemented by mirabait (from the
[MIRA](http://sourceforge.net/projects/mira-assembler/) package). The
names of the reads thus identified, are collected and added to the
list of read names from previous iterations (There is a separate list
for the general baiting and for each thread.) If there are no new
reads identified, then the program stops the iteration, either for the
given thread (specific) or for all the threads (general).

By using a general baiting step before the specific baiting it is
possible to reduce the required run time, because the large input read
file(s) is/are only screened once per iteration and the specific
baiting is confined to screening the reads that are already identified
to be specific during the general baiting.

#### Collecting reads

The identified reads are collected from the read files using
[Seqtk](https://github.com/lh3/seqtk) into internal read files. The
program identifies reads based on the read names, thus the first word
of the identifier line should be unique for each read if they are
single-end or should be the same for both reads of the pair.

#### Assembly

The program can use two assemblers,
[EDENA](http://www.genomic.ch/edena.php) and
[Velvet](https://www.ebi.ac.uk/~zerbino/velvet/), by default, but
there is a skeleton code to add a new assembler to the source code of
the program. Also it is also possible to write [an external perl
script](#using-external_skeleton) that is used by
GRAbB for the assembly. By default
[EDENA](http://www.genomic.ch/edena.php) is used for assembly, but
using command line options the other assemblers can be selected, as
well.

The fact that [single](#modes)- or
[paired-mode](#modes) is selected is passed on to the
assembler program that assembles the specific reads _de novo_. Also,
it is possible to pass additional arguments, such as overlap size
([EDENA](http://www.genomic.ch/edena.php)), to the assembler program
at the invocation of the main program.

#### Testing completion

There are multiple completion criteria that can be specified for the
program. These can be specified for each thread separately
([multi-diff-mode](#modes) or for all of the threads
at once. Also, it is possible to specify multiple criteria for the
same thread or run. In this case the program stops when any one of
these criteria is met.

1. Exhaustive run

    The first completion criterion is implicit, the program stops if
    there is no new information found. This means that either there
    are no new reads found or that the new assembly is identical to
    the bait used for the current iteration. If no other criterion is
    specified then the run is an exhaustive run, since it iterates
    until it cannot find any new information.

2. Length criteria (total length, longest contig or N50)

    This is the simplest explicit completion criterion. There are
    three options that can be used for this setting: total assembly
    size, the length of the longest contig or the N50 value of the
    assembly. This criterion is tested independently for each of the
    threads in [multi-mode](#modes) or, otherwise, for
    the single thread. As mentioned before, multiple criteria can be
    used in a single run, this also applies for the different size
    criteria. These settings are useful when exploring the vicinity of
    a specified sequence region.

3. Matching homologous sequence

     In this case the specific reference sequence is used to identify
     the homologous region within the assembly. To identify the
     matching region, GRAbB uses Exonerate with settings that ensure
     that the whole reference sequence is aligned to the assembly
     contigs. This makes it possible to match sequences that are
     somewhat dissimilar to the reference and may also contain indels,
     causing gaps in the alignment. The possible results of the
     matching can be divided into two groups: the whole sequence is
     matched—also if there are internal regions that correspond to
     gaps in the other sequence—or it is not. In the first case the
     completion criterion is met, and the matched region is extracted
     from the assembly in the same orientation as the reference
     sequence and saved to an output file. In the latter case there
     are two possibilities. If the matched region is larger than the
     one for the previous iteration, then the thread or run will
     continue, but if the size of the match has not improved than the
     thread or run will stop.

#### Modes

1. [Single-](#single) or [paired-mode](#single)

    This information is passed on to the assembler program.

    + Paired-mode

        It is the default if two read files are specified.

        In paired-mode the read files are tested whether GRAbB can
        identify pairs as
        [expected](#collecting-reads).

    + Single-mode

        It is the default mode if there is only one or more than two
        read files specified.
        If two read files are specified than single-mode can be chosen
        by using the [--single](#single) argument.

2. Multi-mode

    Multi-mode is selected by using [--type
    multi](#type) option at invocation. In multi-mode
    the reference file is split into single-entry FASTA files, these
    are referred to as specific references files. Also for each entry
    a separate thread is created. The individual threads are
    independent from each other, thus multiple regions can be
    assembled in a single run without interference from each
    other. The specific reference file is used as initial specific
    bait file. Also in [exonerate-mode](#modes) the
    specific reference file is used for the homology matching.

3. Multi-diff-mode

    Multi-diff-mode is selected by using [--type multi-diff](#type) option at invocation.
It also belongs to the multi-mode with all its properties. The
difference is that when the specific reference files are created if
there is completion criterion specified in the identification line of
the FASTA entry, than this is added to the completion criteria to be
used for the given thread.

4. Exonerate-mode

    Exonerate-mode is selected by using [--type
    exonerate](#type) option at invocation. In this
    case the specific reference sequence is used to identify the
    homologous region within the assembly. To identify the matching
    region, GRAbB uses Exonerate with settings that ensure that the
    whole reference sequence is aligned to the assembly contigs. This
    makes it possible to match sequences that are somewhat dissimilar
    to the reference and may also contain indels, causing gaps in the
    alignment. The possible results of the matching can be divided
    into two groups: the whole sequence is matched—also if there are
    internal regions that correspond to gaps in the other sequence—or
    it is not. In the first case the completion criterion is met, and
    the matched region is extracted from the assembly in the same
    orientation as the reference sequence and saved to an output
    file. In the latter case there are two possibilities. If the
    matched region is larger than the one for the previous iteration,
    then the thread or run will continue, but if the size of the match
    has not improved than the thread or run will stop.

5. Clean-mode

    + Clean-mode

        It is selected by [--clean](#clean).
        GRAbB will remove some internal files to save disk space. But
        there is no information lost, because all the deleted files
        can be reconstructed using the remaining files.

    + Double clean-mode

        It is selected by [--clean --clean](#clean).
        GRAbB will remove some internal files to save disk space.
        At the end of the run all the output files and folders are
        deleted except for the result files.

--------------------------------
#### Arguments

    GRAbB.pl --ref <reference file> --reads <read file 1> [<read file 2>] --folder <directory> --prefix <prefix> [options]

There are four mandatory arguments: [reference
file](#ref), [read file(s)](#reads),
[output folder](#folder) and [prefix for the log and
results](#prefix)

The order of the arguments is not important.

#### ref

    --ref <ref.fas>

The reference file is a FASTA formatted file that contains one or more
sequences. The sequence IDs have to be unique for each sequence (as
required by mirabait). If the file contains multiple sequences and the
program is run in [multi-mode](#modes) then the
reference file is split into separate reference files that contain
only a single sequence, the handling of these files is discussed in
the segment on the [main
loop](#main-loop). Furthermore, the description lines
may contain specification for the completion criterion to be used for
the given sequence that is used if [multi-diff
mode](#modes) is selected. Because the read selection
is based on exact k-mer (31 bp) matching, the reference sequence does
not have to be highly similar to the target sequence.
					      
#### bait

    --bait <bait.fas>

A separate bait file can be specified besides the reference file, this
file together with the reference file will be used as [first
bait](#creating-the-bait). Useful when using special
criterion for the assembly, such as
[homology](#modes).

#### reads

    --reads	<r1.fastq> [<r2.fastq> ...]

Multiple read files can be specified as input. If two read files are
given, then it is assumed that reads are
[paired](#modes), but in
[single-mode](#modes), reads are considered as single
reads. The program identifies read pairs based on the read names, thus
the first word of the identifier line should be the same for both
sequences. The read files may be in FASTA or FASTQ format and may be
compressed (using gnuzip). If [EDENA](http://www.genomic.ch/edena.php)
is selected as assembler program, then all the reads should be of the
same length.

#### folder

    --folder <folder_name>

The directory where all the output will be saved. If the directory is
non-empty then the files it contains can be used like internal
files. In this manner previous runs can be continued, make sure to
remove or replace files that would suggest completion:

Folder structure:

**\<_folder_\>**/

- reference.fas
- bait.fas
- extra_bait.fas
- *\<prefix\>*\_assembly\_thread\_*\<int\>*.fas
- *\<prefix\>*\_result\_thread\_*\<int\>*.fas (if exonerate-mode is selected and the sequence was matched)
- *\<prefix\>*.log
- old\_collection.list
- reads*\<int\>*.fastq or reads*\<int\>*.fasta
- __Round*\<int\>*__/

    + hashstat.bin
    + mirabait.log
    + new\_collection.list
    + positive\_*\<int\>*.txt
    + readpool*\<int\>*.fastq
    + reads*\<int\>*.fastq
    + __thread\_*\<int\>*__/

        * assembly.fas
        * new\_collection.list
        * readpool*\<int\>*.fastq
        * (files or __folders__ generated by the assembler)
        * exonerate.log (if exonerate-mode is selected)
        * result.fas (if exonerate-mode is selected and the sequence was matched)

- __thread\_*\<int\>*__/

    + assembly.fas
    + assembly\_*\<int\>*.fas
    + bait.fas
    + final\_assembly.fas
    + old\_collection.list
    + reference.fas
    + reference.fas.exonerate (if exonerate-mode is selected)
    + result.fas (if exonerate-mode is selected and the sequence was matched)

#### prefix

    --prefix <prefix_of_output>

The prefix for the output files:

* log file
* assembly file
* result file

#### single

    --single

Treat reads as unpaired reads even if two read files are specified

#### min_length

    --min_length=<int>

Minimum size required for a contig to be included for completion
testing and baiting

#### type

    --type <run_type_string>

* When should the iteration stop?

    + The intrinsic criteria are:

        - There are no new reads
      	- There is no assembly
    	- The bait sequence did not change

    + Extra criteria:

        - Length of the assembly

	            --type total=<int>

    	- Length of the longest contig of the assembly

	            --type longest=<int>

    	- N50 value of the assembly

	            --type N50=<int>

    	- Reference matches to assembly (uses exonerate)

	            --type exonerate

* Or to have parallel runs (threads)

        --type multi


* To get the completion criterion for each thread from the reference file
(These criteria overwrite the globally defined
criteria. Write "exhaustive" in the identifier
line if exhaustive should be used.)

    (Need to be used together with multi)

        --type multi_diff


Multiple options can be used at the same time, but then they have to be
typed as a concatenated string (or using quotation)

    --type <multiexonerate | multi_exonerate | "multi exonerate">

#### arg1

    --arg1 "argument1 argument2 argument3"

Arguments needed for the assembler program used for graph/hash generation

#### arg2

    --arg2 "argument1 argument2 argument3"

Arguments needed for the assembler program used for the assembly

#### assembler

    --assembler assembler_name

Specify the assembler to be used:

* [EDENA](http://www.genomic.ch/edena.php)
* [Velvet](https://www.ebi.ac.uk/~zerbino/velvet/)
* [alternative](#adding-to-the-source-code-of-grabb)
* [external](#using-external_skeleton)

#### clean

    --clean

Remove some internal files to save disk space.

    --clean --clean

If specified twice, then only result files are kept, the rest is deleted

--------------------------------
#### Using custom assembler program

GRAbB may use either [EDENA](http://www.genomic.ch/edena.php) or
[Velvet](https://www.ebi.ac.uk/~zerbino/velvet/) to preform the _de
novo_ assembly, by default. There are two possibilities to add other
assemblers to be used by the program.

#### Adding to the source code of GRAbB

First the absolute path or the command for the assembler has to be
added to the source code of GRAbB

    my $alt_1_cmd =   "";
    my $alt_2_cmd =   "";

These lines should be modified so that the command or absolute path of
the executable of the assembler are written between the quotation
marks. The first one is the command used for the first step of the
assembly, the second is the command used for the second step of the
assembly. If the assembly requires only one step, then use [the other
method](#using-external_skeleton) to add the assembler.

In addition, there is a skeleton code within the source code of GRAbB
marked by `# Assemble (Alternative)` at the end of the lines. Within
this block of code there are six lines that need to be changed to add
the new assembler. These lines have exclamation marks at the end.


`$alt_1_cmd` and `$alt_2_cmd` are the commands defined above. The
four command line strings (`$single_first`, `$single_second`,
`$paired_first` and `$paired_second`) have to be completed with
all necessary arguments that are required by the assembler for the
different assembly steps.

* The name of the assembler, this will be displayed by GRAbB during
the run

        my $name = "";

* The code to execute the first step of the assembly for single-end
  sequences (`@arg1`: the arguments specified at invocation after the
  **--arg1** argument; `@$readpool_ref`: list of the readpool files
  created by GRAbB)

        my $single_first = "$alt_1_cmd @arg1 @$readpool_ref";

* The code to execute the second step of the assembly for single-end
  sequences (`@arg2`: the arguments specified at invocation after the
  **--arg2** argument)

        my $single_second = "$alt_2_cmd @arg2";

* The code to execute the first step of the assembly for paired-end
  sequences (`@arg1`: the arguments specified at invocation after the
  **--arg1** argument; `@$readpool_ref`: list of the readpool files
  created be GRAbB)

        my $paired_first = "$alt_1_cmd @arg1 @$readpool_ref";

* The code to execute the second step of the assembly for paired-end
  sequences (`@arg2`: the arguments specified at invocation after the
  **--arg2** argument)

        my $paired_second = "$alt_2_cmd @arg2";

* The relative path of the assembly file created by the assembler. The
  working directory is identical to where the assembly commands have
  been issued

        my $result = "";

After all these modifications are completed, then GRAbB maybe run with
the newly specified assembler program. Using the proper argument:

    --assembler alternative

#### Using external_skeleton

There is a skeleton script that maybe copied and modified to run the
assembly if **external** is selected as assembler.

    #!/usr/bin/perl -w
    use strict;
    
    # A perl script that allows you to use an assembler that is not
      specified in the GRAbB source code
    
    # Get all the parameters
    my ($reads, $paired, $parameters1, $parameters2, $outfile, $format) = @ARGV;
    # Create an array containing the read files
    my @reads = split / /, $reads;
    my @param1 = split / /, $parameters1;
    my @param2 = split / /, $parameters2;
    
    # If there were no extra arguments passed for the assembler at the invocation of GRAbB
    #   then the value passed to this script is "-" ($param1[0] and/or $param2[0])

    # The program to be used for assembly
    my $assmebler = ""; # ADD
    
    # Run the assembly
    `$assmebler @reads`; # ADD
    
    # Specify the expected result file
    my $result = ""; # ADD
    
    # Copy the result file to the expected position
    `perl -ne 'if (/^>/) {print} else {tr/acgtwmrskybvdh/ACGTWMRSKYBVDH/; print}' $result >$outfile`;

The lines marked with `# ADD` have to be adjusted:

* The absolute path or the command for the assembler

        my $assmebler = ""; # ADD

* The command line to execute the assembly. The list of read files is
  stored in the array `@reads`, the arguments specified after
  **--arg1** are stored in the array `@param1`, the arguments
  specified after **--arg2** are stored in the array `@param1`,
  the fact that the assembly is for a single-end or paired-end library
  is stored in the variable `$paired` (the value of the variable is
  either "*single*" or "*paired*") and the variable `$format` stores the
  file format of the read files (either "*fasta*" or "*fastq*").

    ```
`$assmebler @reads`; # ADD
    ```

* The relative path of the assembly file created by the assembler. The
  working directory is identical to where the assembly commands have
  been issued

        my $result = ""; # ADD

After all these modifications are completed and save to the files
_\<external.pl\>_ (the file name can be anything that does not already
exist in the directory), then GRAbB maybe run with the newly specified
assembler program. Using the proper argument:

    --assembler external <external.pl>


--------------------------------
#### Citation

If you use the `GRAbB.pl` or any of the helper programs, please, site the paper describing the program.

_The paper is not yet published_

--------------------------------
#### Contact

Balázs Brankovics <b.brankovics@cbs.knaw.nl>
