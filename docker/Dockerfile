FROM ubuntu:14.04

# Create user and add to sudoers
#  password: grabb
RUN useradd --create-home --shell /bin/bash grabb
RUN echo `echo "grabb\ngrabb\n" | passwd grabb`
RUN usermod -a -G sudo grabb

# Create bin directory for the user
RUN mkdir /home/grabb/bin

# Install couple of prerequisites
RUN echo "deb http://archive.ubuntu.com/ubuntu trusty main universe" > /etc/apt/sources.list
RUN apt-get update
RUN apt-get install -y exonerate
RUN apt-get install -y velvet
RUN apt-get install -y seqtk

# Copy executable files from local folder to the bin folder
COPY edena /home/grabb/bin/
COPY GRAbB.pl /home/grabb/bin/
COPY mirabait /home/grabb/bin/

# Copy helper programs
COPY prinseq-lite.pl /home/grabb/bin/
COPY fasta_shift /home/grabb/bin/
COPY fastq2fasta /home/grabb/bin/
COPY get_overlaps /home/grabb/bin/
COPY interleaved2pairs /home/grabb/bin/
COPY merge_contigs /home/grabb/bin/
COPY pairwise_alignment_score /home/grabb/bin/
COPY rename_fastq /home/grabb/bin/
COPY reverse_complement /home/grabb/bin/
COPY single2pairs /home/grabb/bin/
COPY uniform_length /home/grabb/bin/


# Copy the manual and the tutorial files
COPY README.md /home/grabb/
COPY Examples.md /home/grabb/
COPY Tutorial.md /home/grabb/


# Change the owner of the files
RUN chmod +x /home/grabb/bin/*
RUN chown grabb /home/grabb/bin/*
RUN chown grabb /home/grabb/*

# Make the home directory the working directory
WORKDIR /home/grabb

# Set the user to be the user
USER grabb

# Add the bin folder to the path
ENV PATH /home/grabb/bin:$PATH