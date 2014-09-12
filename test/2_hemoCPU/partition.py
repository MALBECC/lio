#!/usr/bin/python

import fileinput
import random
import sys

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

    for (i,e) in enumerate(array):
        candidates = filter(lambda i: bins[i].slack() >= e,xrange(0,len(bins)))
        if not candidates:
            candidates = [len(bins)]
            bins.append( Bin(capacity) )

        chosen = min(candidates, key=lambda i: bins[i].slack())
        bins[chosen].add(i,e)

    return bins 

def getbestpartition(array, threads):
    left, right = 0, sum(array)+1

    while right-left > 1:
        capacity = (right + left) / 2
        bestbins = getbinsplit(array, capacity)

        if len(bestbins) <= threads:
            right = capacity
        else:
            left = capacity

    return topartition(map(lambda p: p.indexes(), getbinsplit(array, right)))

if __name__ == '__main__':
    array = [int(line) for line in fileinput.input()]
    threads = 4

    partition = getbestpartition(array, threads)
    for (thread, indexes) in partition.items():
        for index in indexes:
            print "%d %d %d" % (thread, index, array[index])
