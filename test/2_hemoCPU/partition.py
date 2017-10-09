#!/usr/bin/python

import random
import sys
import multiprocessing

def topartition(array):
    return dict(map(lambda i: (i, array[i]), xrange(0, len(array))))

class Bin(object):
    def __init__(self, capacity):
        self.capacity = capacity
        self.elements = []

    def add(self, i, e):
        self.elements.append((i, e))

    def indexes(self):
        return map(lambda v: v[0], self.elements)

    def slack(self):
        return self.capacity - sum(map(lambda v: v[1], self.elements))

def getbinsplit(array, capacity):
    bins = [ Bin(capacity) ]

    for (i,e) in array:
        candidates = filter(lambda u: bins[u].slack() >= e,xrange(0,len(bins)))
        if not candidates:
            candidates = [len(bins)]
            bins.append( Bin(capacity) )

        chosen = min(candidates, key=lambda u: bins[u].slack())
        bins[chosen].add(i,e)

    return bins 

def getbestpartition(array, threads):
    left, right = 0, sum([v for u, v in array])+1

    while right-left > 1:
        capacity = (right + left) / 2
        bestbins = getbinsplit(array, capacity)

        if len(bestbins) <= threads:
            right = capacity
        else:
            left = capacity

    return topartition(map(lambda p: p.indexes(), getbinsplit(array, right)))

if __name__ == '__main__':
    array = [int(line) for line in sys.stdin.readlines()]
    array = zip(xrange(0,len(array)), array)
    array = sorted(array, key=lambda p: p[1], reverse=True)
    threads = multiprocessing.cpu_count()
    if len(sys.argv) > 1:
        threads = int(sys.argv[1])

    print "%d %d" % (multiprocessing.cpu_count() / 2 / threads, threads)
    partition = getbestpartition(array, threads)

    costs = dict(array)
    for (thread, indexes) in partition.items():
        for index in indexes:
            print "%d %d %d" % (thread, index, costs[index])
