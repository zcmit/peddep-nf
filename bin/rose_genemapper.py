#!/usr/bin/env python3
"""Fast drop-in for ROSE_geneMapper.py (ROSE 1.3.x).

Runs ROSE's own ROSE_geneMapper.py, with two slow spots of its ROSE_utils
replaced by code that gives the same results:

* makeTranscriptCollection checks each transcript with geneList.count(), a scan
  of the whole gene list per transcript: ~17 min for a GENCODE annotation, paid
  on every geneMapper call. Here the list is turned into a set first.
* geneMapper looks for every TSS within 50 Mb of each region. LocusCollection
  answers that by visiting each 50-bp bin of the 100 Mb window one by one:
  about 4 million dictionary lookups per region, about an hour for a
  genome-wide H3K27ac table. For searches wider than WIDE_BINS bins a sorted
  index is used instead. It returns the same loci in the same order: by strand
  ('+' first), then by the first bin of the search window the locus falls in,
  then by when the locus was added to the collection.

It also fixes a crash in ROSE_geneMapper.py: a region with no gene within 50 Mb
(e.g. on a contig without annotated genes) stopped the script before any output
was written. Such a region now gets an empty CLOSEST_GENE. Where every region
has a gene, the output is the same as ROSE's.

The ROSE module puts this first on PATH as ROSE_geneMapper.py; it finds ROSE's
script further along PATH and takes the same options.
"""
import bisect
import importlib.util
import os
import sys

WIDE_BINS = 2000

# ROSE_geneMapper.py reads distList, which is unset (or left from the previous
# region) when a region has no gene within 50 Mb
GENELESS_FIXES = [
    ("closestGene = allEnhancerGenes[distList.index(min(distList))]\n",
     "closestGene = allEnhancerGenes[distList.index(min(distList))] if allEnhancerGenes else ''\n"),
    ("closestGene = startDict[allEnhancerGenes[distList.index(min(distList))]]['name']\n",
     "closestGene = startDict[allEnhancerGenes[distList.index(min(distList))]]['name'] if allEnhancerGenes else ''\n"),
    ("        signalWithGenes.append([startDict[closestGene]['name'],",
     "        if closestGene: signalWithGenes.append([startDict[closestGene]['name'],"),
]


def find_rose_script(name):
    me = os.path.realpath(__file__)
    for d in os.environ.get('PATH', '').split(os.pathsep):
        path = os.path.join(d or '.', name)
        if os.path.isfile(path) and os.path.realpath(path) != me:
            return path
    sys.exit('%s: ROSE\'s %s is not on PATH' % (os.path.basename(__file__), name))


def load_rose_utils(script):
    """ROSE_utils as ROSE's script would import it: next to the script, else from sys.path."""
    for d in [os.path.dirname(script)] + sys.path:
        path = os.path.join(d or '.', 'ROSE_utils.py')
        if os.path.isfile(path):
            spec = importlib.util.spec_from_file_location('ROSE_utils', path)
            module = importlib.util.module_from_spec(spec)
            sys.modules['ROSE_utils'] = module
            spec.loader.exec_module(module)
            return module
    sys.exit('ROSE_utils.py not found next to %s or on PYTHONPATH' % script)


def patch(rose):
    collection = rose.LocusCollection
    subset_by_bins = collection._LocusCollection__subsetHelper

    def drops_index(method):
        def wrapper(self, *args):
            self.__dict__.pop('_wide_index', None)
            return method(self, *args)
        return wrapper

    def wide_index(self):
        """Per strand key: loci sorted by first bin, as (first bin, last bin, order added, locus)."""
        index = self.__dict__.get('_wide_index')
        if index is None:
            win = self._LocusCollection__winSize
            rows = {}
            for order, lcs in enumerate(self._LocusCollection__loci):
                senses = ['+', '-'] if lcs.sense() == '.' else [lcs.sense()]
                for s in senses:
                    rows.setdefault(lcs.chr() + s, []).append((lcs.start() // win, lcs.end() // win, order, lcs))
            index = {}
            for key, entries in rows.items():
                entries.sort(key=lambda r: r[0])
                index[key] = ([r[0] for r in entries], entries, max(r[1] - r[0] for r in entries))
            self.__dict__['_wide_index'] = index
        return index

    def subset(self, locus, sense):
        win = self._LocusCollection__winSize
        first, last = locus.start() // win, locus.end() // win
        if last - first < WIDE_BINS:
            return subset_by_bins(self, locus, sense)
        sense = sense.lower()
        if ['sense', 'antisense', 'both'].count(sense) != 1:
            raise ValueError("sense command invalid: '" + sense + "'.")
        if locus.sense() == '.' or sense == 'both':
            senses = ['+', '-']
        elif sense == 'sense':
            senses = [locus.sense()]
        else:
            senses = [s for s in ['+', '-'] if s != locus.sense()]
        index = wide_index(self)
        matches = dict()
        for s in senses:
            if locus.chr() + s not in index:
                continue
            starts, entries, span = index[locus.chr() + s]
            lo = bisect.bisect_left(starts, first - span)
            hi = bisect.bisect_right(starts, last)
            hits = [r for r in entries[lo:hi] if r[1] >= first]
            hits.sort(key=lambda r: (max(r[0], first), r[2]))
            for r in hits:
                matches[r[3]] = None
        return list(matches.keys())

    def makeTranscriptCollection(annotFile, upSearch, downSearch, window=500, geneList=[]):
        # as ROSE_utils.makeTranscriptCollection, with a set for the gene list
        if annotFile.upper().count('REFSEQ') == 1:
            refseqTable, refseqDict = rose.importRefseq(annotFile)
            locusList = []
            ticker = 0
            if len(geneList) == 0:
                geneList = list(refseqDict.keys())
            wanted = set(geneList)
            for line in refseqTable[1:]:
                if line[1] in wanted:
                    if line[3] == '-':
                        locus = rose.Locus(line[2], int(line[4]) - downSearch, int(line[5]) + upSearch, line[3], line[1])
                    else:
                        locus = rose.Locus(line[2], int(line[4]) - upSearch, int(line[5]) + downSearch, line[3], line[1])
                    locusList.append(locus)
                    ticker = ticker + 1
                    if ticker % 1000 == 0:
                        print(ticker)
        return rose.LocusCollection(locusList, window)

    collection._LocusCollection__subsetHelper = subset
    collection._LocusCollection__addLocus = drops_index(collection._LocusCollection__addLocus)
    collection.remove = drops_index(collection.remove)
    rose.makeTranscriptCollection = makeTranscriptCollection
    return subset_by_bins


def main():
    script = find_rose_script('ROSE_geneMapper.py')
    patch(load_rose_utils(script))
    with open(script) as f:
        source = f.read()
    for old, new in GENELESS_FIXES:
        source = source.replace(old, new)
    sys.argv[0] = script
    exec(compile(source, script, 'exec'), {'__name__': '__main__', '__file__': script})


if __name__ == '__main__':
    main()
