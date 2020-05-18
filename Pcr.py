#!/usr/bin/env python

'''
classes for sorted lists and binary searching
'''
__author__ = "David Brawand"
__credits__ = ['David Brawand']
__license__ = "MIT"
__version__ = "4.0.0"
__maintainer__ = "David Brawand"
__email__ = "dbrawand@nhs.net"
__status__ = "Production"

import bisect  # to add to sorted lists
import sys

# coordinate segment with metadate
class Segment(tuple):
    def __getitem__(self, index):
        if isinstance(index, tuple):
            return [self[i] for i in index]
        return super(Segment, self).__getitem__(index)

    def __len__(self):
        return abs(self[2]-self[1])


# list of segments (with sorted insertion)
class Segments(list):
    def __init__(self):
        return

    def sortedInsert(self, t):
        bisect.insort(self, t)
        return

    def index(a, x):
        'Locate the leftmost value exactly equal to x'
        i = bisect.bisect_left(a, x)
        if i != len(a) and a[i] == x:
            return i
        raise ValueError

    def find_lt(a, x):
        'Find rightmost value less than x'
        i = bisect.bisect_left(a, x)
        if i:
            return a[i-1]
        raise ValueError

    def find_gt(a, x):
        'Find leftmost value greater than x'
        i = bisect.bisect_right(a, x)
        if i != len(a):
            return a[i]
        raise ValueError

    def find_le(a, x):
        'Find rightmost value less than or equal to x'
        i = bisect.bisect_right(a, x)
        if i:
            return a[i-1]
        raise ValueError

    def find_all_le(a, x):
        'Find rightmost value less than or equal to x and all same ones'
        i = bisect.bisect_right(a, x)
        if i:  # is not beginning of array
            result = [ a[i-1] ]
            while i-2 >= 0 and result[-1][0] == a[i-2][0] and result[-1][1] == a[i-2][1]:
                i -= 1
                result.append(a[i-1])
            return result
        raise ValueError

    def find_ge(a, x):
        'Find leftmost item greater than or equal to x'
        i = bisect.bisect_left(a, x)
        if i != len(a):
            return a[i]
        raise ValueError

    def find_all_ge(a, x):
        'Find leftmost item greater than or equal to x and all same ones'
        i = bisect.bisect_left(a, x)
        if i != len(a):
            result = [ a[i] ]
            while i+1 < len(a) and result[-1][0] == a[i+1][0] and result[-1][1] == a[i+1][1]:
                i += 1
                result.append(a[i])
            return result
        raise ValueError

# Extension of the the segment class
class Primer(Segment):

    def span(self, other):  # max span
        if self[0] == other[0]:
            return max(self[1:3] + other[1:3]) - min(self[1:3] + other[1:3])
        else:  # different chromosome
            return None

if __name__=="__main__":
    seg = Segments()
    for line in sys.stdin:
        f = line.split()
        Segments.sortedInsert(f)
