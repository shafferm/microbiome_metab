"""ko_13_5 parser to get "rare" otus"""

import argparse
from collections import Counter
import cPickle
import numpy
import bisect
from operator import itemgetter

def make_pickle(ko_table):
    f = open(ko_table, 'U')

    #get counts of all otus
    otus = f.readline().strip().split('\t')[1:-1]
    counts = Counter()

    for line in f:
        line = line.strip().split('\t')
        if line[0] == "metadata_KEGG_Description":
            break
        line = line [1:-1]
        for i, ko in enumerate(line):
            if float(ko) > 0:
                counts[i] += 1
    
    #make pickle out of counter for easy reuse to make different quantiles
    output = open("ko_counts.pkl", 'wb')
    cPickle.dump(otus, output)
    cPickle.dump(counts, output)
    output.close()

def take_quantile(quantile):
    pkl = open("ko_counts.pkl", 'rb')
    otus = cPickle.load(pkl)
    counts = cPickle.load(pkl)
    pkl.close()
    #get bottom percentile of counts
    kos = list(counts.items())
    kos.sort(key=itemgetter(1))
    vals = [ko[1] for ko in kos]
    print bisect.bisect(vals, numpy.percentile(vals, quantile))
    kos = kos[:bisect.bisect(vals, numpy.percentile(vals, quantile))]

    del vals
    kos = [ko[0] for ko in kos]
    new_kos = list()
    for ko in kos:
        new_kos.append(otus[ko])

    new_kos = '\n'.join(new_kos)
    f = open("kos_q" + str(int(quantile)) + ".txt", 'w')
    f.write(new_kos)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('to_run', help="pickle or quant")
    parser.add_argument('-q', '--quant', type = float, help="quantile to include")
    parser.add_argument('-f', '--file', default="ko_13_5_precalculated.tab", 
        help="file for unzipped pared kegg file from picrust")
    args = parser.parse_args()
    if args.to_run == "quant" and args.quant != None:
        take_quantile(args.quant)
    if args.to_run == "pickle":
        make_pickle(args.file)