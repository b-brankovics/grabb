FROM brankovics/prerequisites

# Copy executable files from local folder to the bin folder
# Copy helper programs
# Copy the manual and the tutorial files
# Copy the helper script for extra assemblers
COPY ["README.md", "Examples.md", "Tutorial.md", "LICENSE", "external_skeleton.pl", "GRAbB.pl", "helper_programs/fasta_shift", "helper_programs/fastq2fasta", "helper_programs/get_overlaps", "helper_programs/interleaved2pairs", "helper_programs/merge_contigs", "helper_programs/pairwise_alignment_score", "helper_programs/rename_fastq", "helper_programs/reverse_complement", "helper_programs/single2pairs", "helper_programs/uniform_length", "/home/grabb/bin/"]



# Change the owner of the files
RUN chmod +x /home/grabb/bin/* && \
    chown grabb /home/grabb/bin/* && \
    chown grabb /home/grabb/*

# Make the home directory the working directory
WORKDIR /home/grabb

# Set the user to be the user
USER grabb

# Add the bin folder to the path
ENV PATH /home/grabb/bin:$PATH
