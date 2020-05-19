#!/usr/bin/env python
__doc__='''
##################################################################
## ampliconFilter.py | Primer/Amplicon stats for PCR based TAS  ##
## -- QIAGEN PCR TAS protocol                                   ##
## -- filters himera amplicons                                  ##
## -- masks/nulls/soft/hardclips primer sequences               ##
## -- generates primer usage statistics                         ##
##################################################################
'''

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
import datetime, time
from Pcr import Segment, Segments

def readBedPe(fi,genome_fasta):
    genome = pysam.FastaFile(genome_fasta)
    # create primer and amplicon lists
    # FIELDS
    # chrom start
    F, R = Segments(), Segments()
    with open(fi,'r') as fh:
        for i, line in enumerate(fh):
            # ckip commented lines
            if line.startswith('#'):
                continue
            # get fields
            fwd_chrom, fwd_start, fwd_end, rev_chrom, rev_start, rev_end = line.split()
            # FWD primer (index the leftmost position)
            fwd_seq = genome.fetch(fwd_chrom, int(fwd_start), int(fwd_end))
            fwd_segment = Segment([ fwd_chrom, int(fwd_start), i, {
                'coords': (int(fwd_start), int(fwd_end)),
                'seq': fwd_seq
            }])
            F.sortedInsert(fwd_segment)
            # REV primer (index the rightmost position)
            rev_seq = genome.fetch(rev_chrom, int(rev_start), int(rev_end))
            rev_segment = Segment([ rev_chrom, int(rev_end), i, {
                'coords': (int(rev_start), int(rev_end)),
                'seq': revcomp(rev_seq)
            }])
            R.sortedInsert(rev_segment)
    genome.close()
    return F, R

def printClipped(ofh,alnread,primers,globalstat,extratrim,mask,clipping=0,maskadaptor=True):
    '''set quality to zero for leading and trailing primer sequences'''
    # get primer positions and sequence
    newqual, newseq, lprimer, rprimer = '','','',''
    leftprimerstart, leftprimerend = None, None
    riteprimerstart, riteprimerend = None, None

    for b in alnread.aligned_pairs:
        if b[0] is not None and b[1] is not None:
            if primers[0][3]['coords'][0] <= b[1] and b[1] < primers[0][3]['coords'][1]:
                if leftprimerstart is None:
                    leftprimerstart = b[0]
                leftprimerend = b[0]+1
                lprimer += alnread.seq[b[0]]
            elif primers[1][3]['coords'][0] <= b[1] and b[1] < primers[1][3]['coords'][1]:
                if riteprimerstart is None:
                    riteprimerstart = b[0]
                riteprimerend = b[0]+1
                rprimer += alnread.seq[b[0]]

    # adjust trimming or set primer to zero
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
        # change position
        newpos = alnread.pos + leftprimerend
        # add start clipping
        if leftprimerstart != leftprimerend:
            newcigartuples.append((clipping+3,leftprimerend))
        # resolve cigar tuples
        readPos = [0,0]
        for t in alnread.cigartuples:
            if t[0] in (0,1,4,7,8):  # MIS=X
                readPos[0] = readPos[1]
                readPos[1] += t[1]
            else:
                newcigartuples.append((t[0],t[1]))  # no length in read -> preserve
                continue
            # edit/discard/keep cigar tuple
            if readPos[0] < leftprimerend and riteprimerstart < readPos[1]:
                newcigartuples.append((t[0], riteprimerstart-leftprimerend))  # spans both primers -> strip
            elif readPos[1] <= leftprimerend or riteprimerstart <= readPos[0]:
                pass  # cigar within primer -> discard
            elif readPos[0] < leftprimerend and leftprimerend < readPos[1]:
                newcigartuples.append((t[0], readPos[1]-leftprimerend))  # split left
            elif readPos[0] < riteprimerstart and riteprimerstart < readPos[1]:
                newcigartuples.append((t[0], riteprimerstart-readPos[0]))  # split right
            else:
                newcigartuples.append((t[0],t[1]))  # internal -> keep
        # add end clipping
        if riteprimerstart != riteprimerend:
            newcigartuples.append((clipping+3,len(alnread.seq)-riteprimerstart))
        # validate
        try:
            assert sum([t[1] for t in newcigartuples]) == sum([ t[1] for t in alnread.cigartuples ])
        except:
            raise Exception('ClippingError')
        else:
            if clipping > 1:
                newseq = newseq[leftprimerend:riteprimerstart]
                newqual = newqual[leftprimerend:riteprimerstart]
    else:
        newcigartuples = alnread.cigartuples
        newpos = alnread.pos

    # mark if adapter found
    adapterfound = [True if leftprimerstart != 0 else False, True if riteprimerend != len(alnread.seq) else False]

    # collect stats
    loc = 'three' if alnread.is_reverse else 'five'
    globalstat[loc]['found'][int(bool(leftprimerend))] += 1
    globalstat[loc]['mismatch'][matchstring(lprimer, primers[0][3]['seq'][-len(lprimer):]).count('-')] += 1
    globalstat[loc]['missing'][len(primers[0][3]['seq'])-len(lprimer)] += 1
    globalstat[loc]['adapter'][leftprimerstart] += 1

    loc = 'five' if alnread.is_reverse else 'three'
    globalstat[loc]['found'][int(bool(len(alnread.seq)-riteprimerstart))] += 1
    globalstat[loc]['mismatch'][matchstring(rprimer, revcomp(primers[1][3]['seq'])[:len(rprimer)]).count('-')] += 1
    globalstat[loc]['missing'][len(primers[1][3]['seq'])-len(rprimer)] += 1
    globalstat[loc]['adapter'][len(alnread.seq)-riteprimerend] += 1

    # some cheap paranoia code
    try:
        if clipping > 1:
            assert len(newseq) == len(newqual) and \
                len(newseq) == riteprimerstart - leftprimerend
        else:
            assert len(newseq) == len(alnread.seq) and \
                len(newqual) == len(alnread.qual)
    except:
        print(alnread.seq, file=sys.stderr)
        print(newseq, file=sys.stderr)
        print(alnread.qual, file=sys.stderr)
        print(newqual, file=sys.stderr)
        print(alnread.cigar, file=sys.stderr)
        print(newcigartuples, file=sys.stderr)
        print(alnread.pos, file=sys.stderr)
        print(newpos, file=sys.stderr)
        raise

    # update and print
    alnread.seq = newseq
    alnread.qual = newqual
    alnread.cigartuples = newcigartuples
    alnread.pos = newpos
    # add primer names to tag
    try:
        alnread.set_tag('xp', f'{primers[0][2]}|{primers[1][2]}', replace=True)
    except:
        pass
    # write read
    try:
        ofh.write(alnread)
    except:
        print(alnread.qname, alnread.cigarstring, len(alnread.seq), len(alnread.qual), file=sys.stderr)
        print(lprimer, leftprimerstart, leftprimerend, file=sys.stderr)
        print(rprimer, riteprimerstart, riteprimerend, file=sys.stderr)
        raise
    return

'''returns match string with *:match -:mismatch'''
matchstring = lambda x,y: ''.join(['*' if n == y[i] else "-" for i,n in enumerate(x) ])

'''reverse complement'''
revcomp = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A'}[B] for B in x][::-1])

'''is "properly" paired (correct reference and orientation)'''
def is_paired(aln):
    return aln.is_paired and aln.is_reverse != aln.mate_is_reverse
    # return aln.is_proper_pair and aln.is_reverse != aln.mate_is_reverse

'''get matching primers from a list'''
def get_matched(le,ri):
    for l in le:
        for r in ri:
            if l[2] == r[2]:
                return l, r
    return None, None

'''main'''
if __name__=="__main__":
    # timestamp
    ts = time.time()
    # mandatory arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('designfile', help='file with the amplicon design (from Qiagen)')
    parser.add_argument('-m','--metrics', metavar='FILE', help='Metrics output', type=str)
    # I/O
    parser.add_argument('-i','--inputfile', metavar='FILE', help='inputfile (Default: SAM to STDIN)', type=str)
    parser.add_argument('-o','--outputfile', metavar='FILE', help='outputfile (Default: SAM to STDOUT)', type=str)
    parser.add_argument('-d','--discarded', metavar='FILE', help='discarded fragment output', type=str)
    parser.add_argument('-g','--genome', metavar='FILE', help='The indexed genome FASTA file', type=str)

    # other options
    parser.add_argument('--super', help='Allows all primer combinations', action="store_true", default=False)
    parser.add_argument('--mask', help='Sequence hardmasking (Default: False)', action="store_true", default=False)
    parser.add_argument('--clipping', help='Primer clipping [0] noclip, [1] softclip, [2] hardclip (Default: 0)', type=int, default=0)
    parser.add_argument('--primerdistance', metavar='INT', help='maximum distance from nearest possible primer [0]', type=int, default=0)
    parser.add_argument('--tolerance', metavar='INT', help='detection end tolerance [0]', type=int, default=0)
    parser.add_argument('--maxbuffer', metavar='INT', help='maximum read buffer size [100000]', type=int, default=100000)
    parser.add_argument('--extratrim', metavar='INT', help='additional trimming after primer [0]', type=int, default=0)

    args = parser.parse_args()

    # define sequence/quality masking
    if args.mask:
        primermask = [['N','N','N','N'],['!','!','!','!']]
    else:
        primermask = [[],['!','!','!','!']]

    # read design file and create primer libraries (Forward and Reverse)
    print("Reading design file...", end=' ', file=sys.stderr)
    fwd, rev = readBedPe(args.designfile,args.genome)
    print("\rRead design from %s" % args.designfile, file=sys.stderr)

    # open input
    if args.inputfile:
        if args.inputfile[args.inputfile.rfind('.'):] == '.bam':
            infile = pysam.Samfile(args.inputfile,'rb')
        elif args.inputfile[args.inputfile.rfind('.'):] == '.sam':
            infile = pysam.Samfile(args.inputfile,'r')
        else:
            raise Exception('UnkownInputFormat')
    else:
        infile = pysam.Samfile("-",'r')

    # open output
    if args.outputfile:
        if args.outputfile[args.outputfile.rfind('.'):] == '.bam':
            outfile = pysam.Samfile(args.outputfile, 'wb', template = infile )
        elif args.outputfile[args.outputfile.rfind('.'):] == '.sam':
            outfile = pysam.Samfile(args.outputfile, 'w', template = infile )
        else:
            raise Exception('UnkownOutputFormat')
    else:
        outfile = pysam.Samfile("-", 'wh', template = infile )

    # open output
    if args.discarded:
        if args.discarded[args.discarded.rfind('.'):] == '.bam':
            discfile = pysam.Samfile(args.discarded, 'wb', template = infile )
        elif args.discarded[args.discarded.rfind('.'):] == '.sam':
            discfile = pysam.Samfile(args.discarded, 'w', template = infile )
        else:
            raise Exception('UnkownOutputFormat')
    else:
        discfile = None


    # find primer positions
    from collections import defaultdict
    fragment = collections.Counter()
    alnbuffer = {}  # alignment buffer
    gstat = { 'five': collections.defaultdict(collections.Counter), 'three': collections.defaultdict(collections.Counter) }
    for i, aln in enumerate(infile):
        # progressmeter
        if i % 10000 == 0:
            if i % 100000 == 0:
                sys.stderr.write(str(i/1000)+'k')
            else:
                sys.stderr.write('.')
            sys.stderr.flush()
        # filter singletons,unmapped,qcfail,seconday,supplementary
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
        elif aln.is_secondary:
            if discfile:
                try:
                    aln.set_tag('af', 'secondary', replace=True)
                except:
                    pass
                discfile.write(aln)
            pass  # don't print anything secondary
        elif int(aln.flag) & 2048:  #aln.is_supplementary
            if discfile:
                try:
                    aln.set_tag('af', 'supplementary', replace=True)
                except:
                    pass
                discfile.write(aln)
            pass  # don't print anything supplementary
        # buffer and resolve (PairedEnd)
        elif is_paired(aln) and aln.qname in list(alnbuffer.keys()):  # found mate -> print
            # get left right end
            if alnbuffer[aln.qname].is_reverse and not aln.is_reverse:
                firstseg = aln
                lastseg = alnbuffer[aln.qname]
            elif aln.is_reverse and not alnbuffer[aln.qname].is_reverse:
                firstseg = alnbuffer[aln.qname]
                lastseg = aln
            else:
                print(aln.flag, file=sys.stderr)
                raise Exception("ParanoiaGotReal")
            # get ends
            pair_start = firstseg.reference_start + args.tolerance + 1
            pair_end = lastseg.reference_end - args.tolerance - 1
            l = Segment((infile.getrname(firstseg.rname), pair_start))
            r = Segment((infile.getrname(lastseg.rname), pair_end))
            # find primers
            try:
                lefts = fwd.find_all_le(l)
                rites = rev.find_all_ge(r)
                left, rite = (lefts[0], rites[0]) if args.super else get_matched(lefts,rites)
                assert left and rite  # found a matching pair
                assert left[0] == rite[0]  # references match
                assert l[1] - left[3]['coords'][1] < args.primerdistance and rite[3]['coords'][0] - r[1] < args.primerdistance
            except (ValueError, AssertionError, TypeError):  # no primer found or on different chromosomes
                fragment['ectopic'] += 1
                if discfile:
                    try:
                        if left[0] == rite[0]:
                            firstseg.set_tag('af', 'ectopic', replace=True)
                            lastseg.set_tag('af', 'ectopic', replace=True)
                        else:
                            firstseg.set_tag('af', 'primerdistance', replace=True)
                            lastseg.set_tag('af', 'primerdistance', replace=True)
                    except:
                        pass
                    discfile.write(firstseg)
                    discfile.write(lastseg)
            except:
                raise
            else:
                # tag fragment
                if left[2] != rite[2]:
                    # chimera to discard
                    fragment['chimera'] += 1
                    try:
                        firstseg.set_tag('af', 'chimera', replace=True)
                        lastseg.set_tag('af', 'chimera', replace=True)
                    except:
                        pass
                else:
                    # expected amplicon
                    fragment['designed'] += 1

                # write file
                if args.super or left[2] == rite[2]:
                    printClipped(outfile, firstseg, (left, rite), gstat, args.extratrim, primermask, args.clipping)
                    printClipped(outfile, lastseg,  (left, rite), gstat, args.extratrim, primermask, args.clipping)
                else:
                    if discfile:
                        discfile.write(firstseg)
                        discfile.write(lastseg)
            # cleanup buffer
            del alnbuffer[aln.qname]
        elif is_paired(aln):
            # buffer segment
            alnbuffer[aln.qname] = aln

        # not proper pair (both mapped or SingleEnd)
        elif not aln.mate_is_unmapped:  # BOTH MAPPED NOT PROPERLY PAIRED
            if discfile:  # no proper pair flag
                try:
                    aln.set_tag('af', 'notproper', replace=True)
                except:
                    pass
                discfile.write(aln)
        else:
            raise Exception('CantDecideFilter')


        ## check buffer size
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

    # summary statistics
    stats = collections.OrderedDict()

    ## target amplicons
    stats['ANALYSED_INSERTS'] = str(sum(fragment.values()))
    if int(stats['ANALYSED_INSERTS']) > 0:
        stats['FRACTION_ON_TARGET'] = '{:.3f}'.format(fragment['designed']/float(stats['ANALYSED_INSERTS']))
        stats['FRACTION_CHIMERIC'] = '{:.3f}'.format(fragment['chimera']/float(stats['ANALYSED_INSERTS']))
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
