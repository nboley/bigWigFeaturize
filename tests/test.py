import bigWigFeaturize
import pybedtools
import extractors
import os
import cPickle

e = extractors.BigwigExtractor("/users/nasa/lab/dnase-mnase/data/bigwig/GM12878-DNase.bw")
#
#print "Using the bedfile interface"
#n = bigWigFeaturize.new(["/users/nasa/lab/dnase-mnase/data/bigwig/GM12878-DNase.bw"], 20, "/users/nasa/tmp.bed")
#
#print type(n)
#print n
#print n.shape

print "Using the interval list interface"
intervals = [pybedtools.Interval("chr17", 6939053, 6940054)]
i = bigWigFeaturize.new(["/users/nasa/lab/dnase-mnase/data/bigwig/GM12878-DNase.bw"], 1001, intervals=intervals, cache="cache")
print type(i)
print i
print i.shape

#ir = e(intervals)
#print type(ir)
#print ir
#print ir.shape

#import numpy
#numpy.testing.assert_allclose(i, ir)

if os.path.exists("positives.pkl"):
    positives = cPickle.load(open("positives.pkl"))
else:
    positives = []
    for line in open("peaks.bed"):
        parts = line.split()
        positives.append(pybedtools.Interval(parts[0], int(parts[1]), int(parts[2])))
    
    cPickle.dump(positives, open("positives.pkl", "w"))

length = max([i.stop - i.start for i in positives])

print "length", length

pos = bigWigFeaturize.new(["/users/nasa/lab/dnase-mnase/data/bigwig/GM12878-DNase.bw"], length, intervals=positives)
print type(pos)
print pos.shape

pos_cached = bigWigFeaturize.new(["/users/nasa/lab/dnase-mnase/data/bigwig/GM12878-DNase.bw"], length, intervals=positives, cache="cache")
print type(pos_cached)
print pos_cached.shape

posr = e(positives)
print type(posr)
print posr.shape

import numpy
numpy.testing.assert_allclose(pos, posr)
numpy.testing.assert_allclose(pos_cached, posr)

pos2 = bigWigFeaturize.new(["/users/nasa/lab/dnase-mnase/data/bigwig/GM12878-DNase.bw", "/users/nasa/lab/dnase-mnase/data/bigwig/GM12878-DNase.bw"], length, intervals=positives)
print type(pos2)
print pos2.shape
numpy.testing.assert_allclose(pos2, numpy.dstack((posr, posr)))
