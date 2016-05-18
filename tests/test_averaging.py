import bigWigFeaturize
import pybedtools
import extractors
import os
import cPickle
import numpy

e = extractors.BigwigExtractor("/users/nasa/lab/dnase-mnase/data/bigwig/GM12878-DNase.bw")

if os.path.exists("positives.pkl"):
    positives = cPickle.load(open("positives.pkl"))
else:
    positives = []
    for line in open("peaks_average_detail.tsv"):
        parts = line.split()
        positives.append(pybedtools.Interval(parts[0], int(parts[1]), int(parts[2]), score=parts[6]))
    
    cPickle.dump(positives, open("positives.pkl", "w"))

averages = numpy.array([float(i.score) for i in positives]).reshape((len(positives), 1, 1, 1))
print type(averages)
print averages.shape

length = max([i.stop - i.start for i in positives])

print "length", length

pos = bigWigFeaturize.average(["/users/nasa/lab/dnase-mnase/data/bigwig/GM12878-DNase.bw"], intervals=positives)
print type(pos)
print pos.shape

pos_cached = bigWigFeaturize.average(["/users/nasa/lab/dnase-mnase/data/bigwig/GM12878-DNase.bw"], intervals=positives, cache="cache")
print type(pos_cached)
print pos_cached.shape

posr = numpy.mean(e(positives), axis=-1).reshape((len(positives), 1, 1, 1))
print type(posr)
print posr.shape

print "Comparing bigWigAverageOverBed vs naive BigwigExtractor"
numpy.testing.assert_allclose(averages, posr, atol=1e-5)
print "Comparing averaged vs naive BigwigExtractor"
numpy.testing.assert_allclose(pos, posr, atol=1e-5)
print "Comparing averaged cached vs naive BigwigExtractor"
numpy.testing.assert_allclose(pos_cached, posr, atol=1e-5)

pos2 = bigWigFeaturize.average(["/users/nasa/lab/dnase-mnase/data/bigwig/GM12878-DNase.bw", "/users/nasa/lab/dnase-mnase/data/bigwig/GM12878-DNase.bw"], length, intervals=positives)
print type(pos2)
print pos2.shape
stacked = numpy.concatenate((posr, posr), axis=3)
print stacked.shape
numpy.testing.assert_allclose(pos2, stacked, atol=1e-5)
