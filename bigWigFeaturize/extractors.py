import numpy as np
from scipy.io import loadmat

from pybedtools import Interval
from pysam import FastaFile
import wWigIO


class FastaExtractor(object):
    multiprocessing_safe = True

    def __init__(self, datafile, **kwargs):
        self._datafile = datafile

    def __call__(self, intervals, to_mirror=None, **kwargs):
        NUM_SEQ_CHARS = 4
        fasta = FastaFile(self._datafile)

        width = intervals[0].stop - intervals[0].start
        data = np.zeros((len(intervals), 1, NUM_SEQ_CHARS, width))

        for index, interval in enumerate(intervals):
            seq = fasta.fetch(str(interval.chrom), interval.start,
                              interval.stop)
            one_hot_encode_sequence(seq, data[index, 0, :, :])

        # This is performing a reverse complement operation
        if to_mirror is not None:
            for index, mirror in enumerate(to_mirror):
                if mirror:
                    data[index, :, :, :] = data[index, :, ::-1, ::-1]

        return data


def mirror(data, to_mirror):
    for index, mirror in enumerate(to_mirror):
        if mirror:
            data[index, :, :, :] = data[index, :, :, ::-1]


class BigwigExtractor(object):
    multiprocessing_safe = True

    def __init__(self, datafile, **kwargs):
        self._datafile = datafile
        self._halfwidth = kwargs.get('local_norm_halfwidth', None)

    def __call__(self, intervals, to_mirror=None, **kwargs):
        if self._halfwidth:
            width = intervals[0].stop - intervals[0].start
            offset = width // 2 - self._halfwidth

            slopped_intervals = [
                Interval(interval.chrom,
                         interval.start + offset,
                         interval.stop - offset)
                for interval in intervals
            ]

            data = _bigwig_extractor(self._datafile, slopped_intervals,
                                     **kwargs)
            mean = data.mean(axis=3, keepdims=True)
            std = data.std(axis=3, keepdims=True)

            data = (data[:, 0, 0, -offset:-offset + width] - mean) / std
        else:
            data = _bigwig_extractor(self._datafile, intervals, **kwargs)

        if to_mirror is not None:
            mirror(data, to_mirror)

        return data


def _bigwig_extractor(datafile, intervals, **kwargs):
    width = intervals[0].stop - intervals[0].start
    data = np.zeros((len(intervals), 1, 1, width))

    wWigIO.open(datafile)

    for index, interval in enumerate(intervals):
        wWigIO.getData(datafile, interval.chrom, interval.start, interval.stop,
                       data[index, 0, 0, :])

    wWigIO.close(datafile)

    return data


_chr2vplot = {}


class VplotExtractor(object):
    multiprocessing_safe = False

    def __init__(self, datafile, **kwargs):
        self._datafile = datafile
        self._blur_slope = kwargs.get('blur_slope', None)
        self._blur_intercept = kwargs.get('blur_intercept', 1)
        self._blur_order = kwargs.get('blur_order', 0)
        self._blur_precision = kwargs.get('blur_precision', 3)
        self.build_blur_convolutions()
        self._max_fraglen = kwargs.get('max_fraglen', 300)

    def build_blur_convolutions(self, sigma, axis=-1, output=None):
        if self._blur_slope is None:
            self._blur_convolutions = None
            return
        self._blur_convolutions = []
        for i in range(self._max_fraglen):
            sd = float(self._blur_intercept + self._blur_slope*i)
            if sd <= 0:
                continue
            # make the radius of the filter equal to truncate standard deviations
            lw = int(self._blur_precision * sd + 0.5)
            weights = [0.0] * (2 * lw + 1)
            weights[lw] = 1.0
            sum = 1.0
            sd = sd * sd
            # calculate the kernel:
            for ii in range(1, lw + 1):
                tmp = math.exp(-0.5 * float(ii * ii) / sd)
                weights[lw + ii] = tmp
                weights[lw - ii] = tmp
                sum += 2.0 * tmp
            for ii in range(2 * lw + 1):
                weights[ii] /= sum
            # implement first, second and third order derivatives:
            if self._blur_order == 1:  # first derivative
                weights[lw] = 0.0
                for ii in range(1, lw + 1):
                    x = float(ii)
                    tmp = -x / sd * weights[lw + ii]
                    weights[lw + ii] = -tmp
                    weights[lw - ii] = tmp
            elif self._blur_order == 2:  # second derivative
                weights[lw] *= -1.0 / sd
                for ii in range(1, lw + 1):
                    x = float(ii)
                    tmp = (x * x / sd - 1.0) * weights[lw + ii] / sd
                    weights[lw + ii] = tmp
                    weights[lw - ii] = tmp
            elif self._blur_order == 3:  # third derivative
                weights[lw] = 0.0
                sd2 = sd * sd
                for ii in range(1, lw + 1):
                    x = float(ii)
                    tmp = (3.0 - x * x / sd) * x * weights[lw + ii] / sd2
                    weights[lw + ii] = -tmp
                    weights[lw - ii] = tmp
            self._blur_convolutions.append(weights)

    def __call__(self, intervals, to_mirror=None, **kwargs):
        data = _vplot_extractor(self._datafile, intervals, blur_convolutions=self._blur_convolutions, max_fraglen=self._max_fraglen, **kwargs)
        if to_mirror is not None:
            mirror(data, to_mirror)
        return data


def _vplot_extractor(datafile, intervals, blur_convolutions=None, max_fraglen=300, **kwargs):
    # preallocate scratch array for speed
    # assumes all intervals are of the same size - MAY NOT BE TRUE
    vplot_width = intervals[0].stop - intervals[0].start
    data = np.empty((len(intervals), 1, max_fraglen, vplot_width),
                    dtype=np.float32)
    for index, interval in enumerate(intervals):
        try:
            vplot = _chr2vplot[interval.chrom]
        except KeyError:
            matfile = '%s.%s' % (datafile, interval.chrom)
            print 'Loading %s...' % matfile
            vplot = loadmat(matfile)[interval.chrom]
            _chr2vplot[interval.chrom] = vplot

        if blur_convolutions is not None:
            max_blur = max(map(len, blur_convolutions) if len(blur_convolutions) > 0 else 0)
            convolution_count = len(blur_convolutions)
            precision = 3
            scratch = np.empty((max_fraglen, vplot_width + max_blur*precision*2))
            vplot[:max_fraglen, interval.start-precision*max_blur:interval.end+precision*max_blur].todense(
                out=scratch)
            data[index, 0, :, :] = scratch[:, max_blur*precision:vplot_width+max_blur*precision]
            for i, w in blur_convolutions:
                height = max_fraglen - convolution_count + i # blur from the top down (to handle negative intercept)
                convolved = scipy.ndimage.filters.correlate1d(scratch[height], w, -1, 'constant', 0, 0)
                width = len(convolved)
                data[index, 0, height, :] = convolved[width/2 - vplot_width/2:width/2 + vplot_width/2]
        else:
            vplot[:max_fraglen, interval.start:interval.end].todense(
                out=data[index, 0, :, :])
    return data
