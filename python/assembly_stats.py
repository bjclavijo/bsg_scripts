#!/usr/bin/env python3
from statistics import mean, median
from itertools import groupby
from collections import Counter
import sys

def fasta_iter(fasta_name):
    fh = open(fasta_name)

    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))

    for header in faiter:
        # drop the ">"
        headerStr = header.__next__()[1:].strip()

        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())

        yield (headerStr, seq)

class AssemblyStats(object):
    def __init__(self,filename):
        self.lengths_all=[]
        self.lengths_ACGT=[]
        for ff in fasta_iter(filename):
            headerStr, seq = ff
            self.lengths_all.append(len(seq))
            #bplen=0
            #for c in seq:
            #    #if c == 'A' or c == 'C' or c == 'G' or c == 'T' or c == 'a' or c == 'c' or c == 'g' or c == 't':
            #        bplen+=1
            useq=seq.upper()
            bplen=useq.count('A')+useq.count('C')+useq.count('G')+useq.count('T')
            self.lengths_ACGT.append(bplen)
        self.lengths_all.sort(reverse=True)
        self.lengths_ACGT.sort(reverse=True)
    def print_stats(self,points=[0,25,50,75,100]):
        points.sort()
        totalbp=sum(self.lengths_all)
        totalbpACGT=sum(self.lengths_ACGT)
        next_point=0
        current_bp=0
        print ("Total: %s bp"% '{:,}'.format(totalbp))
        for l in self.lengths_all:
            current_bp+=l
            while next_point<len(points) and current_bp>=totalbp*points[next_point]/100:
                print("N%d: %s"%(points[next_point],'{:,}'.format(l)))
                next_point+=1
        print ("Non-ACGT %%: %.2f"%((totalbp-totalbpACGT)*100.0/totalbp))

if __name__ == '__main__':
    for fn in sys.argv[1:]:
        print ("\nStats for %s" % fn)
        AssemblyStats(fn).print_stats()
