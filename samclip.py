#⁄/usr/bin/env python

import re
import sys
import subprocess
import pysam
import argparse
from IO import bamIO
from functools import reduce

"""
Clips primers with samtools primerclip and then reduces clipped primers BQ to avoid realignment/variant calling.

Usage:
    samclip.py primers.bed alignments.bam 

"""

# character to softclipped primers with
QUAL_MASK = '!'
# regex to find query consumin operations in CIGAR string
QUERY_CONSUMERS = re.compile(r'(\d+)([MIS=X])')


def mask_sequence(aln):
    """Masks quality according to changed softclips in cigarstring

    Args:
        aln (pysam.AlignedSegment]): aligned read

    Returns:
        aln (pysam.AlignedSegment]): aligned read with modified BQ in softclipped primers
    """
    try:
        old_cigar = aln.get_tag('OA').split(',')[3]
    except Exception:
        return aln
    new_cigar = reduce(lambda a,x: a+x[1]*int(x[0]), QUERY_CONSUMERS.findall(aln.cigarstring),'')
    old_cigar = reduce(lambda a,x: a+x[1]*int(x[0]), QUERY_CONSUMERS.findall(old_cigar),'')
    # build new base quality string and mask position which are newly softmasked
    new_qual = ''
    for i,q in enumerate(aln.qual):
        new_qual += QUAL_MASK if old_cigar[i] != new_cigar[i] and new_cigar[i] == 'S' else q
    # assign new quality string
    aln.qual = new_qual
    return aln


def clip_primers(args):
    """iterates over all alignments and masks BQ at newly soft-clipped positions

    Args:
        args (ArgmentParser)): parameters as colelcted by argparse
    """
    # build samtools ampliconclip command line
    cmd = ['samtools', 'ampliconclip', '--original', '--soft-clip', '--clipped',
        '--strand', '--no-excluded', '--both-ends']
    if args.rejects:
        cmd += ['--rejects-file', args.rejects,]
    cmd += ['-b', args.bedfile, args.bamfile]

    # run samtools and capture output
    p = subprocess.Popen(' '.join(cmd), stdout=subprocess.PIPE, shell=True)
    alignments = pysam.Samfile(p.stdout, "rb")
    
    # iterare over alignment and write output
    outfile = bamIO(args.outfile, 'w', True, alignments)
    for i, aln in enumerate(alignments):
        # progressmeter
        if i % 10000 == 0:
            if i % 100000 == 0:
                sys.stderr.write(str(i/1000)+'k')
            else:
                sys.stderr.write('.')
            sys.stderr.flush()
        # mask and write
        outfile.write(mask_sequence(aln))


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('bedfile', metavar='BED', help='BED file of primer locations (with strand)')
    parser.add_argument('bamfile', metavar='BAM', help='BAM input file')
    parser.add_argument('-o','--outfile', metavar='FILE', help='outputfile (Default: SAM to STDOUT)', type=str)
    parser.add_argument('-r','--rejects', metavar='FILE', help='Rejected reads (optional)', type=str)
    args = parser.parse_args()

    clip_primers(args)