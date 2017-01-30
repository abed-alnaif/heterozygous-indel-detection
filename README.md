# heterozygous-indel-detection
Matlab script to subtract a reference sequence from a sample sequence -- useful for detecting heterozygous indels.

Description
===========

For sequence traces of heterozygous parts of the genome, some positions
will exhibit "mixed peaks" in which two nucleotides are dominant, since
the two chromosomes differ in sequence at that position. It's easy enough
to detect these locations manually, by examining the sequence trace.
However, this becomes challenging if one of the chromosomes contains an
insertion/deletion (indel), in which case all positions beyond the
location of the indel will exhibit mixed peaks, making it very difficult
to subtract the reference sequence by manual examination.

In cases where the subtracted sequence is ambiguous, IUPAC conventions
are used to indicate the level of ambiguity.

The program uses basepairs 61-80 to align the two sequences. Thus, the 
two sequences should be at least ~100 basepairs long.

USAGE
=====

To use this function, first convert sequence files from to .fcs. If they
are currently .ab1 and the user is using OS X, they can use the following
instructions (modified from this website:
https://www.biostars.org/p/622/; specifically, the method in the post by Jan Van Haarst):

1. Open terminal

2. cd to progs

3. Run convert_trace command; for example: ./convert_trace -out_format scf < salr_EagI_PCR_rev-PAA-04.ab1 > salr_EagI_PCR_rev-PAA-04.fcs

4. move output file from progs directory to desired directory

This function doesn't have any input or output arguments. Simply type the
reference and sample filenames (with the full path, if the files are in a
different directory as this source code) below in the first lines, and 
run the program. The program prints the amount that the sample
sequence needed to be shifted in order to align it with the reference
sequence, and then it prints the sample sequence after the reference
sequence has been subtracted from it. To detect the presence of indels,
one can simply align this output to the reference sequence (using tools
such as BLAST).

Two example .fcs files are provided:

- reference file: 'salr_EagI_PCR_rev-PAA-04.fcs'

- sample file: 'U290_salr_EagI-EagI_for.fcs'

IMPLEMENTATION
==============

First, the two sequences are aligned using basepairs 61-80 of the sample
sequences (we didn't want to use an earlier part of the sequence since
sequencing tends to be less reliable at the start of the target). Then,
the reference sequence is subtracted from the sample sequence.

The program was tested in Matlab R2012b.
