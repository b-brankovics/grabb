Examples of special usages
----------------------
- [Adding Velvet with new paramters to the source code of GRAbB](#adding-velvet-with-new-paramters-to-the-source-code-of-grabb)
- [Using ABySS assembler as external assembler](#using-abyss-assembler-as-external-assembler)
- [Overwriting a previous run directory](#overwriting-a-previous-run-directory)
- [Using a previous run directory](#using-a-previous-run-directory)



----------------------
#### Adding Velvet with new paramters to the source code of GRAbB

This example demonstrates how a new assembler can be added to the source code of `GRAbB.pl`

The following is the work flow used with the assembler that needs to be added to `GRAbB.pl`

    velveth new_velvet 51 -fasta -shortPaired -separate reads1.fas reads2.fas
    velvetg new_velvet -cov_cutoff 50

- The output of this work flow is located at **new_velvet/contigs.fa**
- We want to be able to adjust the kmer used by Velvet at invocation using **--arg1** argument
- We want to be able to adjust the arguments for the velvetg at invocation using **--arg2** argument

Steps: (according to the instruction in [README.md](README.md#adding-to-the-source-code-of-grabb))

1. Add the commands of the assembler

        my $alt_1_cmd = "velveth";
        my $alt_2_cmd = "velvetg";

2. Change the name of the assembler

        my $name = "Velvet adjustable";

3. Specify the command for the single-end graph generation

        my $single_first = "$alt_1_cmd new_velvet @arg1 -$$format_ref -short @$readpool_ref";

4. Specify the command for the single-end assembly

        my $single_second = "$alt_2_cmd new_velvet @arg2";

5.  Specify the command for the paired-end graph generation

        my $paired_first = "$alt_1_cmd new_velvet @arg1 -$$format_ref -shortPaired -separate @$readpool_ref";

6. Specify the command for the paired-end assembly

        my $paired_second = "$alt_2_cmd new_velvet @arg2";

7. Specify the output of the assembler

        my $result = "new_velvet/contigs.fa";

8. Save the modified source code of `GRAbB.pl` to file

Preform the above work flow within `GRAbB.pl`

    ./GRAbB.pl --ref <ref.fas> --folder <folder> --prefix <prefix> --reads <reads1.fas reads2.fas> --assembler alternative --arg1 51 --arg2 "-cov_cutoff 51"


----------------------
#### Using ABySS assembler as external assembler

This example demonstrates how a new assembler can be used with the `external_scaffold.pl` script

The following is the work flow used with the assembler that we want to use

    abyss-pe k=64 name=abyss in='reads1.fasta reads2.fasta'

- The output of this work flow is located at **abyss-contigs.fa**

Steps: (according to the instruction in [README.md](README.md#using-external_scaffold))

1. Add the commands of the assembler

        my $assmebler = "abyss-pe"; # ADD

2. Specify the command for the assembly

        `$assmebler k=64 name=abyss in='@reads'`; # ADD

3. Specify the output of the assembler

        my $result = "abyss-contigs.fa"; # ADD

4. Save the modified source code of external_scaffold.pl to file (e.g. `external_abyss.pl`)

Preform the above work flow within `GRAbB.pl`

    ./GRAbB.pl --ref <ref.fas> --folder <folder> --prefix <prefix> --reads <reads1.fas reads2.fas> --assembler external external_abyss.pl

----------------------
#### Overwriting a previous run directory

How to use the folder (**old**) that has been used for a previous run, but is not needed anymore.

    ./GRAbB.pl --ref data/results/F11_mt.fas --reads data/reads/PH-1_r1.fastq.gz data/reads/PH-1_r2.fastq.gz --folder old --prefix PH-1_mt

`GRAbB.pl` detects that the folder is not empty and asks a few questions, below is an example with the answers needed to continue

    Specified folder (old) is not empty.
    Do you wish to continue?[y/n] y
    Do you want the program to use the content of the folder?
    (Already existing files will not be overwritten, this includes reference, bait and read files. If you wish to replace them delete those files, then run the command again) [y/n] n
    Do you wish to delete the content of the folder?[y/n] y
    Are you sure?[y/n] y

If the previous run was aborted before completion, then `GRAbB.pl` may fail to deleted straightaway. In this case run it again and answer the questions same as before.


----------------------
#### Using a previous run directory

There are two ways to use the folder of a previous run:

1. Rerunning the last Round either starting with baiting or comparing the list of the identified reads of the latest round and the previous one

    1. Rerunning last Round starting with the baiting step

        Some of the files have to be removed:

            rm old/reference.fas old/bait.fas old/extra_bait.fas old/thread_1/reference.fas old/thread_1/bait.fas old/Round#/positive_*.txt old/Round#/readpool*.fast[aq] old/old_collection.list old/thread_1/old_collection.list old/Round#/new_collection.list
        
        - both from the working directory and the thread folders: reference.fas, bait.fas, extra_bait.fas and old_collection.list
        - From the given Round directory: positive_*.txt, readpool*.fast[aq] and new_collcetion.list

        Now `GRAbB.pl` maybe run using the folder and tell `GRAbB.pl` to use the content of the folder

     2. Rerunning last Round starting with the comparing the list of the identified   reads of the latest round and the previous one

       	Some of	the files have to be removed:
       	
       	    rm old/reference.fas old/bait.fas old/extra_bait.fas old/thread_1/reference.fas old/thread_1/bait.fas old/old_collection.list old/thread_1/old_collection.list

        - both from the working directory and the thread folders: reference.fas, bait.fas, extra_bait.fas and old_collection.list

        Now `GRAbB.pl` maybe run using the folder and tell `GRAbB.pl` to use the content of the folder

2. Rerunning `GRAbB.pl` by starting a new round

    Create the next round directory inside the working directory

        mkdir old/Round#

    Create the bait files from the previous assemblies

        for i in old/thread_*; do cp $i/assembly.fas $i/bait.fas; done
        cat old/thread_*/assembly.fas | perl -ne 'if(/\>/){$n++; print ">new$n\n";} else {print}' >old/bait.fas

    Now `GRAbB.pl` maybe run using the folder and tell `GRAbB.pl` to use the content of the folder