#!/usr/bin/env python

'''
Function for input/output
'''
__author__ = "David Brawand"
__credits__ = ['David Brawand']
__license__ = "MIT"
__version__ = "4.0.0"
__maintainer__ = "David Brawand"
__email__ = "dbrawand@nhs.net"
__status__ = "Production"

import pysam

def bamIO(name, mode, stream=False, template=None):
    """
    Inputs:
        file name (optional)
        read or write mode
        get read/write stream if no file name given
        SAM/BAM file template

    BAM/SAM I/O generates appropriate pysam read or write handle (DRY)

    Returns:
        pysam file handle or None
    """
    if name:
        extension = name[name.rfind('.'):]
        if extension == '.bam':
            return pysam.Samfile(name, '{}b'.format(mode), template=template)
        elif extension == '.sam':
            return pysam.Samfile(name, '{}'.format(mode), template=template)
        else:
            raise Exception('UnkownInputFormat')
    else:
        return pysam.Samfile("-",'{}'.format(mode), template=template) if stream else None
