# GFA: Graphical Fragment Assembly (GFA) mode

GFA format is being adapted by assembler programs as one of the
standards for assembly graphs. Assembly graphs can be highly useful
when working with selective assembly.


## GFA format
The purpose of the GFA format is to capture sequence graphs as the
product of an assembly, a representation of variation in genomes,
splice graphs in genes, or even overlap between reads from long-read
sequencing technology.

The GFA format is a tab-delimited text format for describing a set of
sequences and their overlap. The first field of the line identifies
the type of the line. Header lines start with `H`. Segment lines start
with `S`. Link lines start with `L`. A containment line starts with `C`. A
path line starts with `P`.

## Scripts and steps


### `gfa-statistics.pl`

Prints statistics on the **GFA** graphs

### `gfa2fasta.pl`

Converts segments of **GFA** into **FASTA** format file.

### `find-hits.pl`

Runs exonerate and prints segment ids that have a hit.


### `gfa_extract.pl`

Extracts graphs from **GFA** that have segments exactly matching the
ones specified at invocation.

### Steps

1. Convert **GFA** to **FASTA**
2. Find hits
3. Extract graphs
4. Create bait (**FASTA**)
