#!/usr/bin/env python

"""
Clips primers with samtools primerclip and then reduces
clipped primers' BQ to avoid realignment/variant calling.

Usage:
    samclip.py primers.bed alignments.bam 

"""

import re
import sys
import subprocess
import argparse
import pysam
from IO import bamIO

# character to softclipped primers with
QUAL_MASK = '!'
# regex to find query consumin operations in CIGAR string
QUERY_CONSUMERS = re.compile(r'(\d+)([MIS=X])')

def mask_sequence(aln, rstart, rend):
    """Masks quality according to changed softclips in cigarstring

    Args:
        aln (pysam.AlignedSegment): aligned read
        rstart (int): start of non-primer sequence
        rend (int): end of non-primer sequence

    Returns:
        aln (pysam.AlignedSegment): aligned read with modified BQ in softclipped primers
    """
    # get query CIGAR operations
    query_ops = QUERY_CONSUMERS.findall(aln.cigarstring)

    # determine if either end is trimmed by the expected primer pair
    mask_left = rstart == aln.reference_start
    mask_right = aln.reference_end == rend

    # mask BQ if is an expected soft-clipped primer
    if mask_left and query_ops[0][1] == 'S':
        mask_len = int(query_ops[0][0])
        aln.qual = QUAL_MASK*mask_len + aln.qual[mask_len:]
    if mask_right and query_ops[-1][1] == 'S':
        mask_len = int(query_ops[-1][0])
        aln.qual = aln.qual[:-mask_len] + QUAL_MASK*mask_len

    return aln


def resolve(alns,out):
    """resolvers pirmer quality masking for a set of alignments of a read pair

    Args:
        alns ([pysam.AlignedSegments]) List of alignments from a read pair
        out (FileHandle): Writeable file handle for output
    """
    # determine primer boundaries (outermost aligned base after clipping)
    reference_start, reference_end = None, None
    for aln in alns:
        if not reference_start or aln.reference_start < reference_start:
            reference_start = aln.reference_start
        if not reference_end or aln.reference_end > reference_end:
            reference_end = aln.reference_end
    # mask newly softclipped bases in all alignedSegments outside of those boundaries
    for aln in alns:
        # mask and write
        out.write(mask_sequence(aln, reference_start, reference_end))


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
    cmd += ['-b', args.bedfile, args.bamfile, '|', 'samtools', 'sort', '-n']

    # run samtools and capture output
    sys.stderr.write('Running ampliconclip...\n')
    p = subprocess.Popen(' '.join(cmd), stdout=subprocess.PIPE, shell=True)
    alignments = pysam.AlignmentFile(p.stdout, "rb")

    # iterare over alignment and write output
    outfile = bamIO(args.outfile, 'w', True, alignments)
    buffer = []
    for i, aln in enumerate(alignments):
        # progressmeter
        if i % 10000 == 0:
            if i % 100000 == 0:
                sys.stderr.write(str(i/1000)+'k')
            else:
                sys.stderr.write('.')
            sys.stderr.flush()

        # group alignments by read pair
        if not buffer or buffer[-1].qname == aln.qname:
            buffer.append(aln)
        else:
            # resolve
            resolve(buffer, outfile)
            # empty buffer
            buffer = [] 
    # resolve last reads
    resolve(buffer, outfile)


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('bedfile', metavar='BED', help='BED file of primer locations (with strand)')
    parser.add_argument('bamfile', metavar='BAM', help='BAM input file')
    parser.add_argument('-o','--outfile', metavar='FILE', help='outputfile (STDOUT)')
    parser.add_argument('-r','--rejects', metavar='FILE', help='Rejected reads (optional)')
    args = parser.parse_args()

    clip_primers(args)