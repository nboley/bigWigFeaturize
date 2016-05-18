import bigWigFeaturize
import pybedtools
import time
start_time = time.time()
positives = []
for line in open("/srv/scratch/jisraeli/TF_Binding/Sequence_Only_Classification/single_tf_models/single_celltype_models/TCF21/TCF21_IDR_optimal.narrowPeak"):
    parts = line.split()
    positives.append(pybedtools.Interval(parts[0], int(parts[1]), int(parts[2])))

negatives = []
for line in open("/srv/scratch/jisraeli/TF_Binding/Sequence_Only_Classification/single_tf_models/single_celltype_models/TCF21/ATAC_HCASMC_serum_R1.trim.PE2SE.nodup.nonchrM.tn5.pf_summits.bed.500window.notIntersecting0.5TCF21_IDR_optimal"):
    parts = line.split()
    negatives.append(pybedtools.Interval(parts[0], int(parts[1]), int(parts[2])))

print "Interval loading time:", time.time() - start_time

positive_time = time.time()
print "Positives"
n = bigWigFeaturize.new(["/srv/scratch/jisraeli/TF_Binding/Sequence_Only_Classification/single_tf_models/single_celltype_models/TCF21/ATAC_HCASMC_serum_R1.trim.PE2SE.nodup.nonchrM.tn5.pf.fc.signal.bigwig"], 1000, intervals=positives, cache="cache")
print n.shape
print "Positive processing time:", time.time() - positive_time

negative_time = time.time()
print "Negatives"
n = bigWigFeaturize.new(["/srv/scratch/jisraeli/TF_Binding/Sequence_Only_Classification/single_tf_models/single_celltype_models/TCF21/ATAC_HCASMC_serum_R1.trim.PE2SE.nodup.nonchrM.tn5.pf.fc.signal.bigwig"], 340, intervals=negatives, cache="cache")
print n.shape
print "Negative processing time:", time.time() - negative_time
