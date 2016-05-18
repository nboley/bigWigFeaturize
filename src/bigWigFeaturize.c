/* bigWigFeaturize */
/* Copyright (C) 2015 Nasa Sinnott-Armstrong
 * See README in this or parent directory for licensing information. */

#include <stdio.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <Python.h>
#include <numpy/arrayobject.h>

#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "localmem.h"
#include "options.h"
#include "verbose.h"
#include "basicBed.h"
#include "twoBit.h"
#include "bigWig.h"
#include "bits.h"

#define MAX_CACHE_LENGTH 1048

int bedCmpChrom(const void *va, const void *vb)
/* Compare strings such as chromosome names that may have embedded numbers,
 * so that chr4 comes before chr14 */
{
const struct bed *a = *((struct bed **)va);
const struct bed *b = *((struct bed **)vb);
return cmpStringsWithEmbeddedNumbers(a->chrom, b->chrom);
}

struct bed *nextChromInList(struct bed *bedList)
/* Return first bed in list that starts with another chromosome, or NULL if none. */
{
char *chrom = bedList->chrom;
struct bed *bed;
for (bed = bedList->next; bed != NULL; bed = bed->next)
    if (!sameString(bed->chrom, chrom))
        break;
return bed;
}

int simpleFetch(int length, int inputs, struct bbiFile **bbis, struct bed **pBedList, int computeFastaPerBedEntry, float **data, int average)
/* Do the averaging by sorting bedList by chromosome, and then processing each chromosome
 * at once. Faster for long bedLists. */
{
    /* Sort by chromosome. */
    /* It's a bad idea to compute per chromosome but not have things sorted by chromosome. */
    struct bigWigValsOnChrom **chromValsIn = NULL;
    if (!computeFastaPerBedEntry)
    {
        slSort(pBedList, bedCmpChrom);
        chromValsIn = malloc(sizeof(struct bigWigValsOnChrom *) * inputs);
    }
    
    int cell;
    
    long region_length = (length < 0)?-length:length;

    struct bed *bed, *bedList, *nextChrom;

    long region_count = 0;
    for (bedList = *pBedList; bedList != NULL; bedList = bedList->next) region_count++;

    //printf("%d regions with %d signals of length %d\n", region_count, inputs, region_length);

    //printf("allocating...\n");
    *data = malloc(sizeof(float) * inputs * (average?1:region_length) * region_count);
    float *input = *data;
    bzero(input, sizeof(float) * inputs * (average?1:region_length) * region_count);

    struct lm *lm = lmInit(0);

    long successful = 0;
    for (bedList = *pBedList; bedList != NULL; bedList = nextChrom)
    {
        /* Figure out which chromosome we're working on, and the last bed using it. */
        char *chrom = bedList->chrom;
        int fetched = 1;
        if (!computeFastaPerBedEntry)
        {
            verbose(1, "Processing %s\n", chrom);
            nextChrom = nextChromInList(bedList);
            for (cell = 0; cell < inputs; cell++)
                chromValsIn[cell] = bigWigValsOnChromNew();

            for (cell = 0; cell < inputs; cell++)
                fetched = fetched && bigWigValsOnChromFetchData(chromValsIn[cell], chrom, bbis[cell]);
        }
        else
            nextChrom = NULL;
    
        if (fetched)
        {
            /* Loop through beds doing sums and outputting. */
            for (bed = bedList; bed != nextChrom; bed = bed->next)
            {
                chrom = bed->chrom;
                int center = (bed->chromStart + bed->chromEnd)/2;
                int left = center - (region_length/2);
                int right = left + region_length;

                if (average && region_length == 0) {
                    left = bed->chromStart;
                    right = bed->chromEnd;
                }

                if (left < 0) // || right > chromosome_size)
                {
                    fprintf(stderr, "Didn't include %s:%d-%d because the resulting region (%s:%d-%d) overhangs end of chromosome!\n", chrom, bed->chromStart, bed->chromEnd, chrom, left, right);
                    continue;
                }
                if (computeFastaPerBedEntry)
                {
                    for (cell = 0; cell < inputs; cell++)
                    {
                        struct bbiInterval *bigwigdata = bigWigIntervalQuery(bbis[cell], chrom, left, right, lm);
                        struct bbiInterval *current;
                        for (current = bigwigdata; current != NULL; current=current->next) {
                            int start = current->start < left ? left : current->start;
                            int end = current->end > right ? right : current->end;
                            int pos;
                            for (pos = start - left; pos < end - left; pos++)
                                if (average)
                                    input[cell + inputs * successful] += current->val;
                                else
                                    input[pos + region_length * (cell + inputs * successful)] = current->val;
                        }

			if (average)
				input[cell + inputs * successful] /= right-left;
                    }
                }
                else
                {
                    // scale factor is for the averaging
                    double scale = 1.0/(right - left);
                    int pos;
                    for (cell = 0; cell < inputs; cell++)
                        for (pos = left; pos < right; pos++)
                            if (average)
                                input[cell + inputs * successful] += chromValsIn[cell]->valBuf[pos] * scale;
                            else
                                input[pos + region_length * (cell + inputs * successful)] = chromValsIn[cell]->valBuf[pos];
                }
                // only move memory when successful
                successful++;
            }
        } else {
            fprintf(stderr, "No data for chromosome %s!\n", chrom);
            /* If no bigWig data on this chromosome, just output as if coverage is 0 */
            for (bed = bedList; bed != nextChrom; bed = bed->next)
            {
            }
            if (!computeFastaPerBedEntry)
                for (cell = 0; cell < inputs; cell++)
                    bigWigValsOnChromFree(chromValsIn+cell);
        }
    }

    if (!computeFastaPerBedEntry)
        free(chromValsIn);

    return successful;
}

struct bed *getBed(char *regions, PyObject *intervals) {
    if (regions != NULL)
        return bedLoadNAll(regions, 3);
    else {
        int len = PyList_Size(intervals);
        struct bed *bed = malloc(sizeof(struct bed) * len);
        int i;
        for (i = 0; i < len; i++)
        {
            PyObject *entry = PyList_GetItem(intervals, i);
            PyObject *c = PyObject_GetAttr(entry, PyString_FromString("chrom"));
            bed[i].chrom = PyString_AsString(c);
            PyObject *b = PyObject_GetAttr(entry, PyString_FromString("start"));
            bed[i].chromStart = PyInt_AsLong(b);
            PyObject *e = PyObject_GetAttr(entry, PyString_FromString("stop"));
            bed[i].chromEnd = PyInt_AsLong(e);
            bed[i].next = (i < len-1 ? bed+i+1 : NULL) ;
        }
        return bed;
    }
}

PyObject *dataToArray(float *rawdata, int region_count, int signal_count, int region_length, int average)
{
    if (!average) {
    npy_intp *dims = malloc(sizeof(npy_intp) * 4);

    dims[0] = region_count;
    dims[1] = 1;
    dims[2] = signal_count;
    dims[3] = region_length;
  
    PyObject *data = PyArray_SimpleNewFromData(4, dims, NPY_FLOAT32, rawdata);

    Py_INCREF(data);

    return data;
    } else {
    npy_intp *dims = malloc(sizeof(npy_intp) * 4);

    dims[0] = region_count;
    dims[1] = 1;
    dims[2] = 1;
    dims[3] = signal_count;
  
    PyObject *data = PyArray_SimpleNewFromData(4, dims, NPY_FLOAT32, rawdata);

    Py_INCREF(data);

    return data;
    }
}

int cacheFile(char *cache, long hash, int region_length, int average, char *output) {
    return snprintf(output, MAX_CACHE_LENGTH-1, "%s/%lu_%d%s", cache, hash, region_length, average ? "_averaged": "");
}

int saveToCache(float *rawdata, char *cache, long hash, int region_length, int region_count, int signal_count, int average)
{
    char path[MAX_CACHE_LENGTH];
    cacheFile(cache, hash, region_length, average, path);
    //printf("Caching to %s...\n", path);
    // path now contains where it would be were it in cache

    FILE *cachefile = fopen(path, "w");

    if (cachefile == NULL) return -1;

    return fwrite(rawdata, sizeof(float) * (average?1:region_length), region_count * signal_count, cachefile);
}

PyObject *parseArgs(PyObject *signals, char *regions, PyObject *intervals, int length, char *cache, long hash, int average)
{
    struct bed *bedList = getBed(regions, intervals);
    
    int signal_count = PyList_Size(signals);
    
    struct bbiFile **bbis = malloc(sizeof(struct bbiFile *) * signal_count);
    
    int cell;
    for (cell = 0; cell < signal_count; cell++) 
    {
        PyObject *path = PyList_GetItem(signals, cell);
        bbis[cell] = bigWigFileOpen(PyString_AsString(path));
    }
    float *rawdata;
    int successful = simpleFetch(length, signal_count, bbis, &bedList, 1, &rawdata, average);

    if (cache != NULL)
        saveToCache(rawdata, cache, hash, length, successful, signal_count, average);

    return dataToArray(rawdata, successful, signal_count, length, average);
}

long hash_inputs(PyObject *signals, char *regions, PyObject *intervals)
{
    if (regions != NULL)
    {
        return 0;
    }
    else
    {
        srandom(PyList_Size(signals) * PyList_Size(intervals));
        long shash = random();
        int i;

        for (i = 0; i < PyList_Size(signals); i++)
            shash ^= PyObject_Hash(PyList_GetItem(signals, i));

        long hash = random();

        for (i = 0; i < PyList_Size(intervals); i++)
        {
            long ent = PyObject_Hash(PyList_GetItem(intervals, i));
            long h = hash;
            hash ^= ent;
            //printf("Hashed %ld into %ld to get %ld\n", ent, h, hash);
        }

        return hash ^ shash;
    }
}


int inCache(char *cache, long hash, int region_length, int average) {
    char path[MAX_CACHE_LENGTH];
    cacheFile(cache, hash, region_length, average, path);
    // path now contains where it would be were it in cache
    struct stat metadata;
    if (stat(path, &metadata) < 0) {
        return 0;
    } else {
        return S_ISREG(metadata.st_mode);
    }
}

PyObject *getFromCache(char *cache, long hash, int region_length, int region_count, int signal_count, int average)
{
    char path[MAX_CACHE_LENGTH];
    cacheFile(cache, hash, region_length, average, path);
    // path now contains where it would be were it in cache
    struct stat metadata;
    if (stat(path, &metadata) < 0 || !S_ISREG(metadata.st_mode)) {
        perror("stat cache file");
        return NULL;
    }

    // don't trust the bed length; use the actual file size, in case entries were removed
    int successful = metadata.st_size/(sizeof(float) * signal_count * (average?1:region_length));

    int fd = open(path, O_RDONLY);

    if (fd == -1) {
        perror("open cache file");
        return NULL;
    }

    void *rawdata = mmap(0, metadata.st_size, PROT_READ, MAP_SHARED, fd, 0);

    if (rawdata == MAP_FAILED) {
        perror("memmap");
        return NULL;
    }

    return dataToArray(rawdata, successful, signal_count, region_length, average);
}

static PyObject *
bigWigFeaturize_new(self, args, keywds)
    PyObject *self;
    PyObject *args;
    PyObject *keywds;
{  
    PyObject *signals = NULL;
    PyObject *intervals = NULL;
    char *regions = NULL;
    char *cache = NULL;
    int length = 0;
    char average = 0;

    static char *kwlist[] = {"signals", "length", "regions", "intervals", "cache", "average", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "Oi|sOsb", kwlist, 
                                     &signals, &length, &regions, &intervals, &cache, &average))
        return NULL; 

    if ((regions != NULL && intervals != NULL) || (regions == NULL && intervals == NULL)) {
        return NULL;
    }

    long hash = hash_inputs(signals, regions, intervals);
    printf("hashed inputs to %ld\n", hash);
    if (cache != NULL && hash != 0 && inCache(cache, hash, length, average))
    {
        return getFromCache(cache, hash, length, PyList_Size(intervals), PyList_Size(signals), average);
    }
    else
    {
        return parseArgs(signals, regions, intervals, length, cache, hash, average);
    }
}

static PyObject *
bigWigFeaturize_average(self, args, keywds)
    PyObject *self;
    PyObject *args;
    PyObject *keywds;
{  
    PyObject *signals = NULL;
    PyObject *intervals = NULL;
    char *regions = NULL;
    char *cache = NULL;
    int length = 0;

    static char *kwlist[] = {"signals", "length", "regions", "intervals", "cache", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|isOs", kwlist, 
                                     &signals, &length, &regions, &intervals, &cache))
        return NULL; 

    if ((regions != NULL && intervals != NULL) || (regions == NULL && intervals == NULL)) {
        return NULL;
    }

    long hash = hash_inputs(signals, regions, intervals);
    //printf("hashed inputs to %ld\n", hash);
    if (cache != NULL && hash != 0 && inCache(cache, hash, length, 1))
    {
        return getFromCache(cache, hash, length, PyList_Size(intervals), PyList_Size(signals), 1);
    }
    else
    {
        return parseArgs(signals, regions, intervals, length, cache, hash, 1);
    }
}

static PyMethodDef bigWigFeaturize_methods[] = {
    /* The cast of the function is necessary since PyCFunction values
     * only take two PyObject* parameters, and bigWigFeaturize_new() takes
     * three.
     */
    {"new", (PyCFunction)bigWigFeaturize_new, METH_VARARGS|METH_KEYWORDS},
    {"average", (PyCFunction)bigWigFeaturize_average, METH_VARARGS|METH_KEYWORDS},
    {NULL,  NULL}   /* sentinel */
};

PyMODINIT_FUNC
initbigWigFeaturize()
{
  /* Create the module and add the functions */
  Py_InitModule("bigWigFeaturize", bigWigFeaturize_methods);
  import_array();
}
