# Ampliconfilter 
Script to mask primer sequences in aligned DNA amplicons.

Accepts BEDPE files for primer locations and a BAM/SAM file/stream.

Supports filtering of superamplicons/chimeric PCR products.

`ampliconfilter.py --help`  for command line options.

NB: This program implements a standard bisection algorithm to find matching primers.
