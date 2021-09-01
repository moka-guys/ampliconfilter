# Ampliconfilter 
Script to mask primer sequences in aligned DNA amplicons.

Accepts BEDPE files for primer locations and a BAM/SAM file/stream.

Supports filtering of superamplicons/chimeric PCR products.

`ampliconfilter.py --help`  for command line options.

NB: This program implements a standard bisection algorithm to find matching primers.

## Releases

### v1.0

- initial release

### v1.0.1

- fixed an offset error that occured in about 1/1000 reads when the primer sequence contained an indel or was part softclipped
