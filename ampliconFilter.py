#!/usr/bin/env python
__doc__="""
##################################################################
## ampliconFilter.py | Primer/Amplicon stats for PCR based TAS  ##
## -- Paired Primer PCR protocol                                ##
## -- filters chimera amplicons                                 ##
## -- masks/nulls/soft/hardclips primer sequences               ##
## -- generates primer usage statistics                         ##
##################################################################
"""

__author__ = "David Brawand"
__credits__ = ['David Brawand']
__license__ = "MIT"
__version__ = "4.0.0"
__maintainer__ = "David Brawand"
__email__ = "dbrawand@nhs.net"
__status__ = "Production"

import sys
import pysam
import argparse
import collections
import datetime
import time
from collections import defaultdict
from IO import bamIO
from Segment import Segment, Segments
from functools import reduce

def readBedPe(fi,genome_fasta):
    """
    Inputs:
        BEDPE file
        Reference Genome (FASTA) file

    Uses reference genome to obtain primer sequences and create sortable primer segments
    Each Segment object stores primer coordinates and sequences.

    Returns:
        Two Segment objects (segment lists) for forward and reverse primers
    """
    genome = pysam.FastaFile(genome_fasta)
    F, R = Segments(), Segments()
    with open(fi,'r') as fh:
        for i, line in enumerate(fh):
            # skip commented lines
            if line.startswith('#'):
                continue
            # get fields
            fwd_chrom, fwd_start, fwd_end, rev_chrom, rev_start, rev_end = line.split()
            # get FWD primer sequence and create Segment (index the leftmost position)
            fwd_seq = genome.fetch(fwd_chrom, int(fwd_start), int(fwd_end))
            fwd_segment = Segment([ fwd_chrom, int(fwd_start), i, {
                'coords': (int(fwd_start), int(fwd_end)),
                'seq': fwd_seq
            }])
            # Insert created segment to Segments (into sorted coordinate position)
            F.sortedInsert(fwd_segment)
            # get REV primer sequence and create Segment (index the rightmost position)
            rev_seq = genome.fetch(rev_chrom, int(rev_start), int(rev_end))
            rev_segment = Segment([ rev_chrom, int(rev_end), i, {
                'coords': (int(rev_start), int(rev_end)),
                'seq': revcomp(rev_seq)
            }])
            # Insert created segment to Segments (into sorted coordinate position)
            R.sortedInsert(rev_segment)
    genome.close()
    return F, R

def primerClip(alnread,primers,globalstat,extratrim,mask,clipping=0,maskadaptor=True):
    """
    Inputs:
        pysam aligned_pair
        matched primers (list of Segment objects)
        dictionary to collect statistics
        number of bases to trim into biological sequence after primer
        list of list with characters to use for hard masking
        clipping mode (0: no clipping, 1: soft clipping, 2: hard clipping)
        switch to mask adapter sequences (preceding the primer)

    Clips the found primers from a sequence read (hard of softmasking) and writes to file handle
    This will not adjust the mate positions (use sync_mate_pos)

    Returns:
        the primer clipped pysam.AlignedSegment, or None if no sequence remains
    """
    # get primer positions and sequence for read pair
    newqual, newseq, lprimer, rprimer = '','','',''
    leftprimerstart, leftprimerend = None, None
    riteprimerstart, riteprimerend = None, None
    for aligned_pair in alnread.aligned_pairs:  # aligned_pairs [(position in read, position in reference)]
        if aligned_pair[0] is not None and aligned_pair[1] is not None:
            # check if FWD primer overlaps/bookend upstream read boundary
            if primers[0][3]['coords'][0] <= aligned_pair[1] and aligned_pair[1] < primers[0][3]['coords'][1]:
                if leftprimerstart is None:
                    leftprimerstart = aligned_pair[0]
                leftprimerend = aligned_pair[0]+1
                lprimer += alnread.seq[aligned_pair[0]]
            # check if REV primer overlaps/bookend downstream read boundary
            elif primers[1][3]['coords'][0] <= aligned_pair[1] and aligned_pair[1] < primers[1][3]['coords'][1]:
                if riteprimerstart is None:
                    riteprimerstart = aligned_pair[0]
                riteprimerend = aligned_pair[0]+1
                rprimer += alnread.seq[aligned_pair[0]]

    # adjust trimming or set primer start (5') to beginning or end of read
    if leftprimerstart is None:
        leftprimerstart = leftprimerend = 0
    elif extratrim:
        leftprimerend += extratrim

    if riteprimerstart is None:
        riteprimerstart = riteprimerend = len(alnread.seq)
    elif extratrim:
        riteprimerstart -= extratrim

    # build masked sequence
    if mask[0]:
        newseq  = mask[0][0] * leftprimerstart if maskadaptor else alnread.seq[:leftprimerstart]
        newseq += mask[0][1] * (leftprimerend-leftprimerstart)
        newseq += alnread.seq[leftprimerend:riteprimerstart]
        newseq += mask[0][2] * (riteprimerend-riteprimerstart)
        newseq += mask[0][3] * (len(alnread.seq)-riteprimerend) if maskadaptor else alnread.seq[riteprimerend:]
    else:
        newseq = alnread.seq

    # build masked quality string
    if mask[1]:
        newqual  = mask[1][0] * leftprimerstart if maskadaptor else alnread.qual[:leftprimerstart]
        newqual += mask[1][1] * (leftprimerend-leftprimerstart)
        newqual += alnread.qual[leftprimerend:riteprimerstart]
        newqual += mask[1][2] * (riteprimerend-riteprimerstart)
        newqual += mask[1][3] * (len(alnread.qual)-riteprimerend) if maskadaptor else alnread.qual[riteprimerend:]
    else:
        newqual = alnread.qual

    # SOFT/HARDCLIPPING OF PRIMER SEQUENCES
    if clipping:
        newcigartuples = []
        # add start clipping
        if leftprimerstart != leftprimerend:
            # clipping==1 (softclip) => 4 (CSOFT_CLIP) / clipping==2 (hardclip) => 5 (CHARD_CLIP)
            newcigartuples.append((clipping+3, leftprimerend)) 
        # resolve cigar tuples
        readPos = [0,0]
        for t in alnread.cigartuples:
            if t[0] in (0,1,4,7,8):  # MIS=X (query consumers)
                readPos[0] = readPos[1]
                readPos[1] += t[1]
            else:
                newcigartuples.append((t[0],t[1]))  # no length in read -> preserve
                continue
            # edit/discard/keep cigar tuple
            # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*
            # Eg. sequence is split left
            #      |==========================|
            #      readPos[0]        readPos[1]
            #
            #    |===================|------------|===================|
            #    leftprimerstart                  riteprimerstart
            #            leftprimerend                    riteprimerend
            #
            # will become:
            #                        |========|  (left end of cigartuple is clipped)
            #   
            if readPos[0] < leftprimerend and riteprimerstart < readPos[1]:
                # sequence fragment (cigar op) spans both primers -> strip
                newcigartuples.append((t[0], riteprimerstart-leftprimerend))
            elif readPos[1] <= leftprimerend or riteprimerstart <= readPos[0]:
                # sequence fragment within primer -> discard
                pass  
            elif readPos[0] < leftprimerend and leftprimerend < readPos[1]:
                # sequence fragment is split left (preserve OP, modify length)
                newcigartuples.append((t[0], readPos[1]-leftprimerend))
            elif readPos[0] < riteprimerstart and riteprimerstart < readPos[1]:
                # sequence fragment is split right
                newcigartuples.append((t[0], riteprimerstart-readPos[0]))
            else:
                # internal sequence fragment (not primer) -> keep
                newcigartuples.append((t[0],t[1]))
        # add end clipping
        if riteprimerstart != riteprimerend:
            newcigartuples.append((clipping+3,len(alnread.seq)-riteprimerstart))
        # validate
        try:
            assert sum([t[1] for t in newcigartuples]) == sum([ t[1] for t in alnread.cigartuples ])
        except:
            raise Exception('ClippingError')
        else:
            # replace sequence and quality strings if hardclipping
            if clipping > 1:
                newseq = newseq[leftprimerend:riteprimerstart]
                newqual = newqual[leftprimerend:riteprimerstart]
        # change position (needs to be adjusted if begins with softclip on fwd strand)
        # first operation that consumes reference (CIGAR operations MDN=X)
        first_reference_base = 0
        for o, l in newcigartuples:  # (numeric operation code, length of operation) from CIGAR string
            if o in (0,2,3,7,8):  #  break at first reference consumer CIGAR operation (reference consumers)
                break
            elif o in (0,1,4,7,8):  # only advance if query consumed
                first_reference_base += l
        # extract aligned position corresponding to new first reference base
        try:
            # extract reference position of position in query (https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.get_aligned_pairs)
            newpos = alnread.get_aligned_pairs()[first_reference_base][1] 
            assert newpos
        except (AssertionError, IndexError):
            # contains no sequence ofter clipping of primers
            newpos = None
        except:
            print(alnread.qname)
            print(alnread.get_aligned_pairs())
            print(first_reference_base)
            print(alnread.cigartuples, newcigartuples)
            raise Exception("Alignment pair out of range (BAM input valid?)")
    else:
        newcigartuples = alnread.cigartuples
        newpos = alnread.pos

    # collect stats (forward primer)
    loc = 'three' if alnread.is_reverse else 'five'
    # count found primer
    globalstat[loc]['found'][int(bool(leftprimerend))] += 1
    #  count inexact primer match
    globalstat[loc]['mismatch'][matchstring(lprimer, primers[0][3]['seq'][-len(lprimer):]).count('-')] += 1
    # count missing primer
    globalstat[loc]['missing'][len(primers[0][3]['seq'])-len(lprimer)] += 1

    # collect stats (reverse primer)
    loc = 'five' if alnread.is_reverse else 'three'
    globalstat[loc]['found'][int(bool(len(alnread.seq)-riteprimerstart))] += 1
    globalstat[loc]['mismatch'][matchstring(rprimer, revcomp(primers[1][3]['seq'])[:len(rprimer)]).count('-')] += 1
    globalstat[loc]['missing'][len(primers[1][3]['seq'])-len(rprimer)] += 1

    # sanity check if read edit resulted in sequence and quality strings of
    # same length and that the length of the primer has been removed/masked
    try:
        if clipping > 1:
            assert len(newseq) == len(newqual) and \
                len(newseq) == riteprimerstart - leftprimerend
        else:
            assert len(newseq) == len(alnread.seq) and \
                len(newqual) == len(alnread.qual)
    except:
        print(alnread, file=sys.stderr)
        print(alnread.seq, file=sys.stderr)
        print(newseq, file=sys.stderr)
        print(alnread.qual, file=sys.stderr)
        print(newqual, file=sys.stderr)
        print(alnread.cigar, file=sys.stderr)
        print(newcigartuples, file=sys.stderr)
        print(alnread.pos, file=sys.stderr)
        print(newpos, file=sys.stderr)
        raise


    if newpos:
        # update aligned read sequence, quality, CIGAR and position
        alnread.seq = newseq
        alnread.qual = newqual
        alnread.cigartuples = newcigartuples
        alnread.pos = newpos
        # add primer names to tag
        try:
            alnread.set_tag('xp', f'{primers[0][2]}|{primers[1][2]}', replace=True)
        except:
            pass
        # return primerclipped read
        return alnread
    else:
        return None

"""returns match string with *:match -:mismatch"""
matchstring = lambda x,y: ''.join(['*' if n == y[i] else "-" for i,n in enumerate(x) ])

"""reverse complement"""
revcomp = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A'}[B] for B in x][::-1])

"""flatten list of lists"""
flatten = lambda x: [item for sublist in x for item in sublist]

"""check if all elements in list are the same"""
all_same = lambda x: not x or x.count(x[0]) == len(x)

def sync_mate_positions(alns):
    """
    Inputs:
        a list of pysam.AlignedSegments

    Checks for query name match and adjust mate positions respectively

    Returns:
        Nothing (mutates input)
    """
    if len(set(map(lambda a: a.qname, alns))) == 1:
        for i in range(len(alns)):
            if not alns[i].is_secondary:
                for j in range(len(alns)):
                    if i != j and alns[i].is_reverse != alns[j].is_reverse:
                        alns[j].pnext = alns[i].pos
                        if not alns[j].is_secondary:
                            alns[i].pnext = alns[j].pos

    return alns

def is_paired(aln):
    """
    Input:
        pysam.AlignedSegment

    "properly" paired (correct reference and orientation)
    
    NB: This is used instead of the properly_paired flag as rare superamplicons
    fall outside the expected insert size distribution and are not marked as such. 
    
    NB: It does not check if the orientation is FR as RF orientations are discarded due
    to no biological sequence.

    Returns:
        Boolean
    
    """
    return aln.is_paired and aln.is_reverse != aln.mate_is_reverse

def is_complete(alns):
    """
    Input:
        list of pysam.AlignedSegment
    
    returns true if all secondary and supplementary alignments are included (if any) and
    forward and reverse primary alignments exists

    Returns:
        Boolean

    """
    try:
        # check if at least two
        assert len(alns)>1
        # check if forward and reverse primary alignment present
        rev_primary = [aln.is_reverse for aln in alns if not aln.is_secondary]
        assert any(rev_primary) and not all(rev_primary)
        # check if seconday and if all there
        secondary_alignments = flatten([ aln.get_tag('SA').rstrip(';').split(';') \
            for aln in alns if not aln.is_secondary and aln.has_tag('SA')])
        assert len(secondary_alignments) == len([aln for aln in alns if aln.is_secondary])
    except AssertionError:
        return False
    return True

def get_matched_primer_pair(le,ri):
    """
    Inputs:
        two lists of primers

    get matching primers from a list (get primer pair from same amplicon)
    
    Returns:
        matched primer pair (same amplicon) if any

    """
    for l in le:
        for r in ri:
            if l[2] == r[2]:
                return l, r
    return None, None

def primersOverlap(l,r):
    """
    Inputs:
        left primer Segment
        right primer Segment

    Check if two primers overlap

    Returns:
        Boolean
    """
    return l[3]['coords'][1] > r[3]['coords'][0]

def expectedPrimerPair(l,r):
    """
    Inputs:
        left primer Segment
        right primer Segment

    Check if 2 primers come from the same primer pair

    Returns:
        Boolean
    """
    return l[2] == r[2]

if __name__=="__main__":
    """
    Main CLI interface
    """
    # timestamp
    ts = time.time()
    # mandatory arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('designfile', metavar='BEDPE', help='BEDPE file with the amplicon design')
    parser.add_argument('-m','--metrics', metavar='FILE', help='Metrics output', type=str)
    # I/O
    parser.add_argument('-i','--inputfile', metavar='FILE', help='inputfile (Default: SAM to STDIN)', type=str)
    parser.add_argument('-o','--outputfile', metavar='FILE', help='outputfile (Default: SAM to STDOUT)', type=str)
    parser.add_argument('-d','--discarded', metavar='FILE', help='discarded fragment output', type=str)
    parser.add_argument('-g','--genome', metavar='FILE', help='The indexed genome FASTA file', type=str)

    # other options
    parser.add_argument('--super', help='Allows all primer combinations (chimera/superamplicons)', action="store_true", default=False)
    parser.add_argument('--mask', help='Sequence hardmasking (Default: False)', action="store_true", default=False)
    parser.add_argument('--clipping', help='Primer clipping [0] noclip, [1] softclip, [2] hardclip (Default: 0)', type=int, choices=[0,1,2], default=0)
    parser.add_argument('--primerdistance', metavar='INT', help='maximum distance from nearest possible primer 5\' end [0]', type=int, default=0)
    parser.add_argument('--tolerance', metavar='INT', help='detection end tolerance [0]', type=int, default=0)
    parser.add_argument('--maxbuffer', metavar='INT', help='maximum read buffer size [100000]', type=int, default=100000)
    parser.add_argument('--extratrim', metavar='INT', help='additional trimming after primer [0]', type=int, default=0)

    args = parser.parse_args()

    # define sequence/quality masking characters (for sequence and quality string)
    # [ fwd_adapter, fwd_primer, rev_primer, rev_adapter ]
    if args.mask:
        primermask = [['N','N','N','N'],['!','!','!','!']]
    else:
        # we mask the quality in any case as this prevents variantcallers from using primer sequences through realignment of softclipped sequences
        primermask = [[],['!','!','!','!']]

    # read design file and create primer libraries (Forward and Reverse)
    print("Reading design file...", end=' ', file=sys.stderr)
    fwd, rev = readBedPe(args.designfile,args.genome)
    print("\rRead design from %s" % args.designfile, file=sys.stderr)

    # open SAM/BAM input file/stream
    infile = bamIO(args.inputfile, 'r', True)

    # open output for filtered/clipped reads
    outfile = bamIO(args.outputfile, 'w', True, infile)

    # open output stream for discarded reads
    discfile = bamIO(args.discarded, 'w', False, infile)

    # fragment counter and primer counters
    fragment = collections.Counter()
    gstat = {'five': collections.defaultdict(collections.Counter),
             'three': collections.defaultdict(collections.Counter)}

    # alignment buffer (to find and combine read pairs)
    alnbuffer = defaultdict(list)  # alignment buffer

    # find primer positions
    for i, aln in enumerate(infile):
        # progressmeter
        if i % 10000 == 0:
            if i % 100000 == 0:
                sys.stderr.write(str(i/1000)+'k')
            else:
                sys.stderr.write('.')
            sys.stderr.flush()
        # filter singletons,unmapped,qcfail,seconday,supplementary (write to discarded file)
        if not aln.is_unmapped and aln.mate_is_unmapped:  # M-
            if discfile:
                try:
                    aln.set_tag('af', 'singleton', replace=True)
                except:
                    pass
                discfile.write(aln)
            pass  # don't print any singletons
        elif aln.is_unmapped:  # -M or --
            if discfile:
                try:
                    aln.set_tag('af', 'unmapped', replace=True)
                except:
                    pass
                discfile.write(aln)
            pass  # don't print anything unmapped
        elif aln.is_qcfail:
            if discfile:
                try:
                    aln.set_tag('af', 'qc', replace=True)
                except:
                    pass
                discfile.write(aln)
            pass  # don't print anything that failed QC
        elif aln.is_supplementary:  # supplementary are non-linear
            if discfile:
                try:
                    aln.set_tag('af', 'supplementary', replace=True)
                except:
                    pass
                discfile.write(aln)
            pass  # don't print anything supplementary
        # read is paired and has a buffered mate
        elif is_paired(aln):
            # buffer segment
            alnbuffer[aln.qname].append(aln)
        elif not aln.mate_is_unmapped:
            # discard not proper pair (both mapped or SingleEnd)
            if discfile:  # no proper pair flag
                try:
                    aln.set_tag('af', 'notproper', replace=True)
                except:
                    pass
                discfile.write(aln)
        else:
            raise Exception('CantDecideFilter')

        # check if read pair (and secondary alignemnts) are all in buffer and can be resolved
        if is_complete(alnbuffer[aln.qname]):
            alns = alnbuffer[aln.qname]

            # check if primary alignments on same chromosome, else discard as non-linear
            primary_same_reference = all_same([ a.reference_name for a in alns])
            if not primary_same_reference:
                # discard all non primary with other reference
                if discfile:  # no proper pair flag
                    for a in alns:
                        try:
                            a.set_tag('af', 'nonlinear', replace=True)
                        except:
                            pass
                        discfile.write(a)
                del alnbuffer[aln.qname]
                continue

            # get leftmost and rightmost alignment (which are relevant for primer matching)
            firstseg, lastseg = None, None
            for a in alns:
                if a.is_reverse:
                    if not lastseg or (lastseg.reference_end < a.reference_end):
                        lastseg = a
                else:
                    if not firstseg or (firstseg.reference_start > a.reference_start):
                        firstseg = a
            # check FR orientation
            try:
                assert firstseg and lastseg
            except Exception as e:
                raise Exception("Reads are not in FR orientation")

            # get ends of read pair (add tolerance offset for finding primer)
            pair_start = firstseg.reference_start + args.tolerance + 1
            pair_end = lastseg.reference_end - args.tolerance - 1
            # create Segments for each read
            l = Segment((infile.getrname(firstseg.rname), pair_start))
            r = Segment((infile.getrname(lastseg.rname), pair_end))
            # find all matching primers by bisection
            try:
                # get all primers that match (primers in each list are all the same but from difference amplicon pairs)
                lefts = fwd.find_all_le(l)
                rites = rev.find_all_ge(r)
                # if superamplicons(chimera) allowed use first primers found
                # else only get primers from matching pairs
                left, rite = (lefts[0], rites[0]) if args.super else get_matched_primer_pair(lefts,rites)
                assert left and rite  # found a matching pair
                assert left[0] == rite[0]  # references match
                # fwd read starts before primer 3'end or up to primerdistance after
                assert l[1] - left[3]['coords'][1] < args.primerdistance
                # same logic for reverse primer
                assert rite[3]['coords'][0] - r[1] < args.primerdistance
            except (ValueError, AssertionError, TypeError):
                # no primer found or on different chromosomes
                # tag alignments with reason and write to discard file
                fragment['ectopic'] += 1
                if discfile:
                    try:
                        if left[0] == rite[0]:
                            for a in alns:
                                a.set_tag('af', 'ectopic', replace=True)
                        else:
                            for a in alns:
                                a.set_tag('af', 'primerdistance', replace=True)
                    except:
                        pass
                    for a in alns:
                        discfile.write(a)
            except:
                raise
            else:
                # get read pair primer info
                expected = expectedPrimerPair(left, rite)
                overlaps = primersOverlap(left, rite)
                # tag fragment and count
                if not expected:
                    # superamplicon (chimera)
                    fragment['chimera'] += 1
                    try:
                        firstseg.set_tag('af', 'chimera', replace=True)
                        lastseg.set_tag('af', 'chimera', replace=True)
                    except:
                        pass
                elif overlaps:
                    # no biological sequence
                    fragment['short'] += 1
                    try:
                        firstseg.set_tag('af', 'short', replace=True)
                        lastseg.set_tag('af', 'short', replace=True)
                    except:
                        pass
                else:
                    # expected amplicon
                    fragment['designed'] += 1

                # write file if no primer overlaps and expected or super
                if not overlaps and (args.super or expectedPrimerPair(left, rite)):
                    # clip primers from read alignments, return None if no sequence remains
                    clipped_alns = list(map(lambda a: primerClip(a, (left, rite), gstat,args.extratrim, primermask, args.clipping), alns))
                    if all(clipped_alns):
                        # synchronise new mate pos to primary alignments
                        clipped_alns = sync_mate_positions(alns)
                        # print all
                        for alnread in clipped_alns:
                            try:
                                outfile.write(alnread)
                            except Exception as e:
                                print(alnread, file=sys.stderr)
                                raise Exception("Could not write clipped read")
                    else:
                        # discard pair if one sequence only contains primer
                        if discfile:
                            firstseg.set_tag('af', 'primeronly', replace=True)
                            lastseg.set_tag('af', 'primeronly', replace=True)
                            discfile.write(firstseg)
                            discfile.write(lastseg)

                else:
                    if discfile:
                        discfile.write(firstseg)
                        discfile.write(lastseg)
            
            # remove read from alignment buffer when resolved
            del alnbuffer[aln.qname]

        # delete empty created buffer entry
        if not len(alnbuffer[aln.qname]):
            del alnbuffer[aln.qname]

        ## check buffer overflow (avoid memory usage escalation)
        if len(alnbuffer) > args.maxbuffer:
            print("ERROR: Buffer overflow (%d)" % args.maxbuffer, file=sys.stderr)
            raise Exception("BufferOverflow")


    # CLOSE OPEN FILEHANDLES
    if discfile:
        discfile.close()
    outfile.close()
    infile.close()
    print('[DONE]', file=sys.stderr)

    # check if all printed as pairs
    try:
        assert len(alnbuffer) == 0
    except:
        print("Could not find all mate pairs! Missed ", len(list(alnbuffer.keys())), file=sys.stderr)
        print("First 50 mate-less mappings:\n", '\n'.join(list(alnbuffer.keys())[:50]), file=sys.stderr)
        #raise Exception("TruncatedInput")

    # calculate summary statistics
    stats = collections.OrderedDict()
    stats['ANALYSED_INSERTS'] = str(sum(fragment.values()))
    if int(stats['ANALYSED_INSERTS']) > 0:
        stats['FRACTION_ON_TARGET'] = '{:.3f}'.format(fragment['designed']/float(stats['ANALYSED_INSERTS']))
        stats['FRACTION_CHIMERIC'] = '{:.3f}'.format(fragment['chimera']/float(stats['ANALYSED_INSERTS']))
        stats['FRACTION_SHORT'] = '{:.3f}'.format(fragment['short']/float(stats['ANALYSED_INSERTS']))
        stats['FRACTION_ECTOPIC'] = '{:.3f}'.format(fragment['ectopic']/float(stats['ANALYSED_INSERTS']))
    else:
        stats['FRACTION_ON_TARGET'] = 'NA'
        stats['FRACTION_CHIMERIC'] = 'NA'
        stats['FRACTION_ECTOPIC'] = 'NA'

    # print the metrics
    outfh = open(args.metrics,'w') if args.metrics else sys.stderr
    print("## custom.ampliconfilter.metrics", file=outfh)
    print("#", ' '.join(sys.argv), file=outfh)
    print("## custom.ampliconfilter.metrics", file=outfh)
    print("# Started on:", datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S'), file=outfh)
    print("\n# FILTER METRICS", file=outfh)
    print('\t'.join(list(stats.keys())), file=outfh)
    print('\t'.join(list(stats.values())), file=outfh)
    print("\n# PRIMER METRICS", file=outfh)
    print('\t'.join(['seqend','stat','len','count']), file=outfh)
    for end in sorted(gstat.keys()):
        for stat in sorted(gstat[end].keys()):
            for length in sorted(gstat[end][stat].keys()):
                print('\t'.join(map(str,[end,stat,length,gstat[end][stat][length]])), file=outfh)
    print(file=outfh)  # add empty line
    if args.metrics:
        outfh.close()
