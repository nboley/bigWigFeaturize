import bigWigFeaturize
import pybedtools
import bigWig
import time
start_time = time.time()
positives = []
for line in open("big.bed"):
    parts = line.split()
    positives.append(pybedtools.Interval(parts[0], int(parts[1]), int(parts[2])))

print "Interval loading time:", time.time() - start_time

positive_time = time.time()
print "Positives"
bw = bigWig.BigwigBinaryExtractor(["/mnt/data/epigenomeRoadmap/signal/consolidated/macs2signal/foldChange/E123-DNase.fc.signal.bigwig"], cache="cache", local_norm_halfwidth=5000)
n = bw(positives)
#bigWigFeaturize.new(["/mnt/data/epigenomeRoadmap/signal/consolidated/macs2signal/foldChange/E123-DNase.fc.signal.bigwig"], 1000, intervals=positives, cache="cache")
print n.shape
print "Positive processing time:", time.time() - positive_time
