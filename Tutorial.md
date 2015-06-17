Tutorials for GRAbB
----------------------
- [Reconstruction of the mitochondrial genome using the mitogenome of a related species](#reconstruction-of-the-mitochondrial-genome-using-the-mitogenome-of-a-related-species)
- [Reconstruction of the mitochondrial genome using multiple references](#reconstruction-of-the-mitochondrial-genome-using-multiple-references)
- [Reconstruction of the rDNA region using partial sequences as references](#reconstruction-of-the-rdna-region-using-partial-sequences-as-references)
- [Extract homologous sequences](#extract-homologous-sequences)
- [Extract homologous sequence with polymorphic allele](#extract-homologous-sequence-with-polymorphic-allele)
- [Length criterion](#length-criterion)
- [Exhaustive assembly of a linear sequence](#exhaustive-assembly-of-a-linear-sequence)
- [Complex run with diverse criteria](#complex-run-with-diverse-criteria)
- [Reconstruction of the rDNA region from actual reads](#reconstruction-of-the-rdna-region-from-actual-reads)



----------------------
#### Reconstruction of the mitochondrial genome using the mitogenome of a related species

- This tutorial shows how to assemble the mitogenome of *Fusarium graminearum* (PH-1) by using *F. oxysporum* (F11) as a reference

Files used

- reads: **data/reads/PH-1_r1.fastq.gz** **data/reads/PH-1_r2.fastq.gz**
- reference: **data/references/F11_mt.fas**

Steps

1. Run `GRAbB.pl` by specifying the above files

        ./GRAbB.pl --ref data/results/F11_mt.fas --reads data/reads/PH-1_r1.fastq.gz data/reads/PH-1_r2.fastq.gz --folder tutorial1 --prefix PH-1_mt

    - The program will run for 9 iterations
    - Ending with an assembly (**tutorial1/PH-1_mt_assembly_thread_1.fas**) with a single contig

2. Run `get_overlaps` script to see if the contig overlaps with itself

        ./helper_programs/get_overlaps tutorial1/PH-1_mt_assembly_thread_1.fas 15

    - 15 specifies the minimum overlap to be detected
    - The script will display which contig overlaps with which contig and what is the size of the overlap
    - The script shows that the contig overlaps with itself and the overlap is 84 bp long

3. Circularize the genome

        ./helper_programs/merge_contigs tutorial1/PH-1_mt_assembly_thread_1.fas 15 circular.fas

    - This cuts the 84 bp overlap from one of the ends of the contig and marks the entry as circular
    - The resulting sequence is saved to circular.fas

4. Visually compare the sequence with the provided result sequence

        nucmer data/references/NC_009493.fas circular.fas; mummerplot out.delta

    - For this `mummer` has to be installed on the system (Otherwise you can use `exonerate` to align the two sequences)
    - With this method we can see if the two sequences are in the same orientation or not or whether they are in the same phase
    - In this case our sequence has to be reverse complemented and shifted

5. Reverse complement the sequence

        ./helper_programs/reverse_complement circular.fas >reverse.fas

6. Shift the sequence into the proper phase

    - For this we either need a short sequence which marks where the sequence should start or we have to specify the position

    a. Find the position where the sequence should start

             exonerate data/references/NC_009493.fas reverse.fas | less

	- Look for the following line "Query range: 0 ->"
	- At the beginning of the alignment we can see which location corresponds to the start position in our sequence (89100)

            ./helper_programs/fasta_shift -i reverse.fas -p 89100 >final.fas

    b. Create a fasta file which contains the start of the PH-1 mitogenome (data/references/NC_009493.fas)

             head -2 data/references/NC_009493.fas >ref.fas
             ./helper_programs/fasta_shift -i reverse.fas -r ref.fas >final.fas

7. Compare result with the correct sequence

        ./helper_programs/pairwise_alignment_score final.fas data/results/PH-1_mt.fas

    - This script outputs a comparison of the two sequences, the percentage marks the similarity of the two
    - We can see that two sequences are 100% identical

----------------	
#### Reconstruction of the mitochondrial genome using multiple references

- This tutorial shows how to assemble the mitogenome using multiple reference sequences
- Also in this case the downstream process is more complicated than in the first tutorial

Files used

- reads: **data/reads/F11_r1.fastq.gz** **data/reads/F11_r2.fastq.gz**
- reference: **data/references/published_mts.fas**

Steps

1. Run `GRAbB.pl` by specifying the above files

        ./GRAbB.pl --ref data/references/published_mts.fas --reads data/reads/F11_r1.fastq.gz data/reads/F11_r2.fastq.gz --folder tutorial2 --prefix F11_mt

    - The program will run for 2 iterations
    - Ending with an assembly (**tutorial2/F11_mt_assembly_thread_1.fas**) with 1 contig

2. Run `get_overlaps` script to see if the contigs overlap

        ./helper_programs/get_overlaps tutorial2/F11_mt_assembly_thread_1.fas 15

    - 15 specifies the minimum overlap to be detected
    - The script will display which contig overlaps with which contig and what is the size of the overlap
    - The script shows that the contig overlaps with itself and the overlap is 24 bp long

3. Circularize the genome

        ./helper_programs/merge_contigs tutorial2/F11_mt_assembly_thread_1.fas 15 circular.fas

    - This cuts the 24 bp overlap from one of the ends of the contig and marks the entry as circular
    - The resulting sequence is saved to circular.fas

4. Visually compare the sequence with the provided result sequence

        nucmer data/results/F11_mt.fas circular.fas; mummerplot out.delta

    - For this `mummer` has to be installed on the system (Otherwise you can use `exonerate` to align the two sequences)
    - With this method we can see if the two sequences are in the same orientation or not or whether they are in the same phase
    - In this case our sequence has to be shifted, but it doesn't have to be reversed

5. Shift the sequence into the proper phase

    - For this we either need short sequence which marks where the sequence should start or we have to specify the position (See also [tutorial 1](#reconstruction-of-the-mitochondrial-genome-using-the-mitogenome-of-a-related-species))
    - Create a fasta file which contains the start of the provided result sequence (**data/results/F11_mt.fas**)

            head -2 data/results/F11_mt.fas >ref.fas
            ./helper_programs/fasta_shift -i circular.fas -r ref.fas >final.fas

6. Compare result with the correct sequence

        ./helper_programs/pairwise_alignment_score final.fas data/results/F11_mt.fas

    - This script outputs a comparison of the two sequences, the percentage marks the similarity of the two
    - We can see that two sequences are 100% identical

---------------
#### Reconstruction of the rDNA region using partial sequences as references

- The rDNA region is not circular DNA, but the sequence is present in multiple direct repeats so it behaves like a circular sequence in the computational sense
- In this tutorial we are going to use IGS and ITS sequences as a refernce to assemble the rDNA region

Files used

    - reads: **data/reads/PH-1_r1.fastq.gz** **data/reads/PH-1_r2.fastq.gz**
    - reference: **data/references/HQ651161.fas**
    - bait: **data/references/JQ423253.fas**

Steps

1. Run `GRAbB.pl` by specifying the above files

        ./GRAbB.pl  --ref data/references/HQ651161.fas --bait data/references/JQ423253.fas --reads data/reads/PH-1_r1.fastq.gz data/reads/PH-1_r2.fastq.gz --folder tutorial3 --prefix rDNA

    - The program will run for 4 iterations
    - Ending with an assembly (**tutorial3/rDNA_assembly_thread_1.fas**) with a single contig

2. Run `get_overlaps` script to see if the contig overlaps with itself

        ./helper_programs/get_overlaps tutorial3/rDNA_assembly_thread_1.fas 15

    - 15 specifies the minimum overlap to be detected
    - The script will display which contig overlaps with which contig and what is the size of the overlap
    - The script shows that the contig overlaps with itself and the overlap is 100 bp long

3. Circularize the genome

        ./helper_programs/merge_contigs tutorial3/rDNA_assembly_thread_1.fas 15 circular.fas

    - This cuts the 100 bp overlap from one of the ends of the contig and marks the entry as circular
    - The resulting sequence is saved to circular.fas

4. Visually compare the sequence with the provided result sequence

        nucmer data/results/PH-1_rDNA.fas circular.fas; mummerplot out.delta

    - For this mummer has to be installed on the system (Otherwise you can use `exonerate` to align the two sequences)
    - With this method we can see if the two sequences are in the same orientation or not or whether they are in the same phase
    - In this case our sequence has to be shifted, but it doesn't have to be reversed

5. Shift the sequence into the proper phase

    - For this we either need short sequence which marks where the sequence should start or we have to specify the position (See also [tutorial 1](#reconstruction-of-the-mitochondrial-genome-using-the-mitogenome-of-a-related-species))

    - Create a fasta file which contains the start of the provided result sequence (**data/results/PH-1_rDNA.fas**)

        head -2 data/results/PH-1_rDNA.fas >ref.fas
        ./helper_programs/fasta_shift -i circular.fas -r ref.fas >final.fas

6. Compare result with the correct sequence

        ./helper_programs/pairwise_alignment_score final.fas data/results/PH-1_rDNA.fas

    - This script outputs a comparison of the two sequences, the percentage marks the similarity of the two
    - We can see that two sequences are 100% identical

------------------	
#### Extract homologous sequences

- In this tutorial three partial gene sequences will be used to extract homologous sequences from another strain
- Reference sequences from Foc1 strain and reads from F11 strain (both them are _F. oxysporum_)

Files used

    - reads: **data/reads/F11\_r1.fastq.gz** **data/reads/F11\_r2.fastq.gz**
    - reference: **data/references/Foc1\_three\_loci.fas**

Steps

1. Run `GRAbB.pl` by specifying the above files

        ./GRAbB.pl --ref data/references/Foc1_three_loci.fas --reads data/reads/F11_r1.fastq.gz data/reads/F11_r2.fastq.gz --folder tutorial4 --prefix F11_three_loci --type multi_exonerate

    - The program runs for one iteration and saves the matching regions into separate files (__tutorial4/F11\_three\_loci\_result\_*\[123\]*.fas__)

2. Compare results with the correct sequences

        ./helper_programs/pairwise_alignment_score tutorial4/F11_three_loci_result_thread_1.fas data/results/F11_cal_360-620.fas
        ./helper_programs/pairwise_alignment_score tutorial4/F11_three_loci_result_thread_2.fas data/results/F11_rpb2_rc1640-2400.fas
        ./helper_programs/pairwise_alignment_score tutorial4/F11_three_loci_result_thread_3.fas data/results/F11_tef1a_679-1239.fas

    - This script outputs a comparison of the two sequences, the percentage marks the similarity of the two
    - We can see that each pair of sequences are 100% identical to each other

-------------------
#### Extract homologous sequence with polymorphic allele

- In this tutorial part of the rDNA region will be extracted from a _Fusarium graminearum_ (PH-1) strain using a reference from a _F. oxysporum_ (F11) strain

    + The selected region has highly similar sequences at both termini whithin the two species while the middle part of the sequences are highly dissimilar

- As an extra bait the whole rDNA sequence of the _F. graminearum_ strain will be specified, thus the program only has to search the specified read files only once
- This method can be used when extracting a region that is highly variable or has multiple possible alleles (both has somewhat conserved termini)

    + By specifying an extra bait file with all known possible alleles the run can be made more efficient, since the input read files have to be searched less often
    + This also means that if we don't know which allele our strain has we can still extract it with any of the known alleles

Files used

    - reads: **data/reads/PH-1\_r1.fastq.gz** **data/reads/PH-1\_r2.fastq.gz**
    - reference: **data/references/F11\_rDNA\_5000-7705.fas**
    - bait: **data/results/PH-1\_rDNA.fas**

Steps

1. Run `GRAbB.pl` by specifying the above files

        ./GRAbB.pl --ref data/references/F11_rDNA_5000-7705.fas --bait data/results/PH-1_rDNA.fas --reads data/reads/PH-1_r1.fastq.gz data/reads/PH-1_r2.fastq.gz --folder tutorial5 --prefix PH-1_IGS --type multi_exonerate

    - The program runs for one iteration and saves the matching region into separate file (**tutorial5/PH-1\_IGS\_result\_thread_1.fas)
    - The program is run without the extra bait file then the program runs for two iterations

2. Compare result with the correct sequence

        ./helper_programs/pairwise_alignment_score tutorial5/PH-1_IGS_result_thread_1.fas data/results/PH-1_rDNA_5001-7733.fas

    - This script outputs a comparison of the two sequences, the percentage marks the similarity of the two
    - We can see that two sequences are 100% identical
    - The correct sequence was identified using blasting the reference sequences against **data/results/PH-1\_rDNA.fas**

-------------
#### Length criterion

- This tutorial demonstrates how to use one of the length criteria (total assembly length) for the assembly
- This can be useful if the vicinity of the sequence is unknown

Files used

    - reads: **data/reads/F11\_r1.fastq.gz** **data/reads/F11\_r2.fastq.gz**
    - reference: **data/references/AY527619.fas**

Steps

1. Run `GRAbB.pl` by specifying the above files

        ./GRAbB.pl --ref data/references/AY527619.fas --reads data/reads/F11_r1.fastq.gz data/reads/F11_r2.fastq.gz --folder tutorial6 --prefix F11_rns_vicintiy --type total=3000

    - The program runs for 3 iterations
    - Ending with an assembly (**tutorial6/F11\_rns\_vicinity\_assembly\_thread\_1.fas**) with a single contig
    - The size of the contigs is 3856 bp

2. Compare result with the correct sequence

        ./helper_programs/pairwise_alignment_score tutorial6/F11_rns_vicintiy_assembly_thread_1.fas data/results/F11_mt_39896-43751.fas

    - This script outputs a comparison of the two sequences, the percentage marks the similarity of the two
    - We can see that two sequences are 100% identical
    - The correct sequence was identified using blasting the result sequences against **data/results/F11\_mt.fas**

---------
#### Exhaustive assembly of a linear sequence

- In this tutorial we are going to use a partial sequence of the tef1a gene to extract the whole gene from the assembly
- In this dataset the tef1a sequence is a linear genomic element that has no continuation, this makes it similar to linear virus genomes

    + tef1a would be further extend to its surroundings in real life sample

Files used

    - reads: **data/reads/Foc1\_r1.fastq.gz** **data/reads/Foc1\_r2.fastq.gz**
    - reference: **data/references/AY527535.fas**

Steps

1. Run `GRAbB.pl` by specifying the above files

        ./GRAbB.pl --ref data/references/AY527535.fas --reads data/reads/Foc1_r1.fastq.gz data/reads/Foc1_r2.fastq.gz --folder tutorial7 --prefix Foc1_tef1a

    - The program will run for 4 iterations
    - Ending with an assembly (**tutorial7/Foc1\_tef1a\_assembly\_thread\_1.fas**) with a single contig

2. Visually compare the sequence with the provided result sequence

        nucmer tutorial7/Foc1_tef1a_assembly_thread_1.fas data/results/Foc1_tef1a.fas; mummerplot out.delta

    - For this `mummer` has to be installed on the system (Otherwise you can use `exonerate` to align the two sequences)
    - With this method we can see if the two sequences are in opposite orientations

        + Thus our sequence has to reverse complemented

3. Reverse complement the sequence

        ./helper_programs/reverse_complement tutorial7/Foc1_tef1a_assembly_thread_1.fas >final.fas

4. Compare result with the correct sequence

        ./helper_programs/pairwise_alignment_score final.fas data/results/Foc1_tef1a.fas

    - This script outputs a comparison of the two sequences, the percentage marks the similarity of the two
    - We can see that two sequences are 100% identical
	
--------------
#### Complex run with diverse criteria

- In this tutorial the different runs that were discussed will be combined into a single run

    + Extract homologous sequence with polymorphic allele (IGS)
    + Extract homologous sequences (cal, rpb2, tef1a)
    + Length criterion (mitochondrial genome total=3000)
    + Exhaustive assembly of a linear sequence
    + Reconstruction of the mitochondrial genome using the mitogenome of a related species

Files used

- reads: **data/reads/Foc1\_r1.fastq.gz** **data/reads/Foc1\_r2.fastq.gz**
- reference: **data/references/tutorial8\_reference.fas**

    + The identifier lines contain the criterion specifications (see brackets after each file name)
    + Contains the sequences from the following files

        * __data/results/PH-1\_rDNA\_5001-7733.fas__     (exonerate)
        * __data/references/Foc1\_cal\_360-620.fas__     (exonerate)
        * __data/results/F11\_rpb2\_rc1640-2400.fas__    (exonerate)
        * __data/references/Foc1\_tef1a_680-1240.fas__  (exonerate)
        * __data/references/AY527619.fas__             (total=3000)
        * __data/references/AY527535.fas__
        * __data/results/PH-1\_mt.fas__

- bait: **data/references/tutorial8\_bait.fas**

    + Contains the sequences from the following files

        * __data/references/F11\_rDNA\_5000-7705.fas__
        * __data/results/PH-1\_rDNA\_5001-7733.fas__

Steps

1. Run `GRAbB.pl` by specifying the above files

        ./GRAbB.pl --ref data/references/tutorial8_reference.fas --bait data/references/tutorial8_bait.fas --reads data/reads/Foc1_r1.fastq.gz data/reads/Foc1_r2.fastq.gz --folder tutorial8 --prefix complex --type multi_diff


    - The program will run for 7 iterations
    - The program reads the completion criterion from each sequence entry in the reference file
    - Ending with 7 assemblies files (**tutorial8/complex\_assembly\_thread\_\[*1234567*\].fas**) and 4 results files (**tutorial8/complex\_result\_thread\_\[*1234*\].fas**)

        + The assemblies contain one contig each

2. Run `get_overlaps` script to see if the mitochondrial contig overlaps with itself

        ./helper_programs/get_overlaps tutorial8/complex_assembly_thread_7.fas 15

    - 15 specifies the minimum overlap to be detected
    - The script will display which contig overlaps with which contig and what is the size of the overlap
    - The script shows that the contig overlaps with itself and the overlap is 24 bp long

3. Circularize the mitochondrial genome

        ./helper_programs/merge_contigs tutorial8/complex_assembly_thread_7.fas 15 circular.fas

    - This cuts the 24 bp overlap from one of the ends of the contig and marks the entry as circular
    - The resulting sequence is saved to circular.fas

4. Visually compare the sequences with the reference sequences

    - For this `mummer` has to be installed on the system (Otherwise you can use `exonerate` to align the two sequences)
    - With this method we can see if the two sequences are in the same orientation or not or whether they are in the same phase

            nucmer data/results/Foc1_mt_38516-42174.fas tutorial8/complex_assembly_thread_5.fas; mummerplot out.delta

    - In this case the sequence is in the same orientation

            nucmer data/results/Foc1_tef1a.fas tutorial8/complex_assembly_thread_6.fas; mummerplot out.delta

    - In this case the sequence has to be reverse complemented

            nucmer data/results/PH-1_mt.fas circular.fas; mummerplot out.delta

    - In this case the sequence has to be reverse complemented and shifted

5. Reverse complement the sequences

        ./helper_programs/reverse_complement tutorial8/complex_assembly_thread_6.fas >reverse_thread6.fas
        ./helper_programs/reverse_complement circular.fas >reverse.fas

6. Shift the sequence into the proper phase

        head -2 data/results/Foc1_mt.fas >ref.fas
        ./helper_programs/fasta_shift -i reverse.fas -r ref.fas >final.fas

    - For this we either need short sequence which marks where the sequence should start or we have to specify the position
    - Create a fasta file which contains the start of the Foc1 mitogenome (**data/results/Foc1\_mt.fas**)
    - Then shift the sequence using the reference file (**ref.fas**)

7. Compare result with the correct sequence

        ./helper_programs/pairwise_alignment_score tutorial8/complex_result_thread_1.fas data/results/Foc1_rDNA_5000-7871.fas
        ./helper_programs/pairwise_alignment_score tutorial8/complex_result_thread_2.fas data/references/Foc1_cal_360-620.fas
        ./helper_programs/pairwise_alignment_score tutorial8/complex_result_thread_3.fas data/references/Foc1_rpb2_rc1640-2400.fas
        ./helper_programs/pairwise_alignment_score tutorial8/complex_result_thread_4.fas data/references/Foc1_tef1a_680-1240.fas
        ./helper_programs/pairwise_alignment_score tutorial8/complex_assembly_thread_5.fas data/results/Foc1_mt_38516-42174.fas      
        ./helper_programs/pairwise_alignment_score reverse_thread6.fas data/results/Foc1_tef1a.fas
        ./helper_programs/pairwise_alignment_score final.fas data/results/Foc1_mt.fas

    - This script outputs a comparison of the two sequences, the percentage marks the similarity of the two
    - We can see that each pair of sequences are 100% identical

------------------	
#### Reconstruction of the rDNA region from actual reads

- In this tutorial the first 100000 read pairs will be used from the SRR550150 sequencing run file

    + This is produced from _Fusarium oxysporum_ f. sp. _cubense_ race 1 strain Foc1

Files used

    - reads: **data/reads/SRR550150\_100000\_1.fastq.gz** **data/reads/SRR550150\_100000\_2.fastq.gz**
    - reference: **data/results/Foc1_rDNA.fas**

Steps

1. Run `GRAbB.pl` by specifying the above files

        ./GRAbB.pl --ref data/results/Foc1_rDNA.fas --reads data/reads/SRR550150_100000_1.fastq.gz data/reads/SRR550150_100000_2.fastq.gz --folder tutorial9 --prefix Foc1_rDNA

    - The program will run for 2 iterations
    - Ending with an assembly (**tutorial9/Foc1\_rDNA\_assembly\_thread\_1.fas**) with a single contig

2. Run `get_overlaps` script to see if the contig overlaps with itself

        ./helper_programs/get_overlaps tutorial9/Foc1_rDNA_assembly_thread_1.fas 15

    - 15 specifies the minimum overlap to be detected
    - The script will display which contig overlaps with which contig and what is the size of the overlap
    - The script shows that the contig overlaps with itself and the overlap is 34 bp long

3. Circularize the genome

        ./helper_programs/merge_contigs tutorial9/Foc1_rDNA_assembly_thread_1.fas 15 circular.fas

    - This cuts the 34 bp overlap from one of the ends of the contig and marks the entry as circular
    - The resulting sequence is saved to circular.fas

4. Visually compare the sequence with the provided result sequence

        nucmer data/results/Foc1_rDNA.fas circular.fas; mummerplot out.delta

    - For this `mummer` has to be installed on the system (Otherwise you can use `exonerate` to align the two sequences)
    - With this method we can see if the two sequences are in the same orientation or not or whether they are in the same phase
    - In this case our sequence has to be shifted, but it doesn't have to be reverse complemented

5. Shift the sequence into the proper phase

        head -2 data/results/Foc1_rDNA.fas >ref.fas
        ./helper_programs/fasta_shift -i circular.fas -r ref.fas >final.fas

    - For this we either need short sequence which marks where the sequence should start or we have to specify the position (See also [tutorial 1](#reconstruction-of-the-mitochondrial-genome-using-the-mitogenome-of-a-related-species))
    - Create a fasta file which contains the start of the provided result sequence (**data/results/Foc1\_rDNA.fas**)

6. Compare result with the correct sequence

        ./helper_programs/pairwise_alignment_score final.fas data/results/Foc1_rDNA.fas

    - This script outputs a comparison of the two sequences, the percentage marks the similarity of the two
    - We can see that two sequences are 100% identical
