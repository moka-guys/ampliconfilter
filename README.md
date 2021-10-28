# Ampliconfilter 
Scripts to mask primer sequences in aligned DNA amplicons.

## ampliconfilter.py

Accepts BEDPE files for primer locations and a BAM/SAM file/stream.

Supports filtering of superamplicons/chimeric PCR products.

`ampliconfilter.py --help`  for command line options.

NB: This program implements a standard bisection algorithm to find matching primers.

## samclip.py (alternative from v1.0.2)

This script wraps the samtools amplicon clipping tool. It will reduce quality of newly soft-clipped bases (primer sequences) to prevent realignment and use in calling variants.

NB: the samtools driven clipping doesnt consider secondary alignments to be part of the same amplicon and therefore does some overclipping

### Usage

`samclip.py -o output.bam -r rejects.bam primers.bed input.bam`

The BED files _must_ contain the intervals of primers and contain strand information.
```
1	1234	1270	amp_1_fwd	.	+
1	1390	1412	amp_1_rev	.	-
```
A script to convert the BEDPE file from ampliconfilter is included `convert_bedpe_to_bed.sh input.bedpe > output.bed`.


# Releases

### v1.0

- initial release

### v1.0.1

- fixed an offset error that occured in about 1/1000 reads when the primer sequence contained an indel or was part softclipped

### v1.0.2

- handling of secondary alignments by ampliconfilter
- added samtools based clipping with quality score reduction on softclipped bases (samclip.py)

