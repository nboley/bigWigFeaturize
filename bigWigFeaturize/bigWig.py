import bigWigFeaturize
import pybedtools

def encode_peaks_bigwig_into_array(peaks, bigwig_fnames, cache=None, average=False):
    lengths = [len(p) for p in peaks]

    if len(set(lengths)) != 1:
        raise ValueError("Can't run with different lengths of sequences! Lengths %s" % set(lengths))
    else:
        length = max(lengths)

    if average:
        if isinstance(cache, basestring):
            return bigWigFeaturize.new(bigwig_fnames, length-1, intervals=peaks, cache=cache, average=average)
        else:
            return bigWigFeaturize.new(bigwig_fnames, length-1, intervals=peaks, average=average)
    else:
        if isinstance(cache, basestring):
            return bigWigFeaturize.new(bigwig_fnames, length-1, intervals=peaks, cache=cache)
        else:
            return bigWigFeaturize.new(bigwig_fnames, length-1, intervals=peaks)

class BigwigBinaryExtractor(object):
    multiprocessing_safe = True

    def __init__(self, datafiles, **kwargs):
        self._datafiles = datafiles
        self._halfwidth = kwargs.get('local_norm_halfwidth', None)
        self._cache = kwargs.get('cache', None)

    def __call__(self, intervals, to_mirror=None, **kwargs):
        if self._halfwidth:
            width = intervals[0].stop - intervals[0].start
            offset = width // 2 - self._halfwidth

            slopped_intervals = [
                pybedtools.Interval(interval.chrom,
                         interval.start + offset,
                         interval.stop - offset)
                for interval in intervals
            ]

            data = encode_peaks_bigwig_into_array(slopped_intervals, self._datafiles, cache=self._cache)
            mean = data.mean(axis=3, keepdims=True)
            std = data.std(axis=3, keepdims=True)

            data = (data[:, 0, 0, -offset:-offset + width] - mean) / std
        else:
            data = encode_peaks_bigwig_into_array(intervals, self._datafiles, cache=self._cache)

        if to_mirror is not None:
            mirror(data, to_mirror)

        return data

def mirror(data, to_mirror):
    for index, mirror in enumerate(to_mirror):
        if mirror:
            data[index, :, :, :] = data[index, :, :, ::-1]
