"""otu_ko_network.py
Takes metagenome_contributions.py output file from PICRUSt and uses it to creat a otu, KO metabolic
network, includes functionality to only include seed sets based on a user defined quantile.
Outputs files which can be loaded directly into Cytoscape for visualization of networks.  Filtering
by pathway or by significance to be added later in addition to including compounds that are used in
reactions by KO's.
"""

import argparse
from collections import Counter
import numpy
import bisect
from operator import itemgetter
from network import *
from parse_KEGG import KEGG_Parser
from make_gg_genomes import parse_gg

kegg_parser = KEGG_Parser()

def parse_meta_contribs(contribs_loc):
    """
        input: contribs_loc = location of metagenome contributions file
        outputs:    genomes = dict with OTU's as keys and sets of ko's as values
    """
    f = open(contribs_loc, 'U')
    f.readline() #toss header
    genomes = dict()
    for line in f:
        line = line.strip().split('\t')
        if line[2] in genomes:
            genomes[line[2]].add(line[0])
        else:
            genomes[line[2]]=set([line[0]])
    f.close()
    return genomes

def filter_genes(genes_to_keep, genomes):
    """takes a list of genes to keep and removes all others from the genome, takes in a set of
    genes so the input gene data structure must be cast to a set.
    """
    new_genomes = dict()
    for genome in genomes:
        genes = genes_to_keep & genomes[genome]
        if len(genes) > 0:
            new_genomes[genome] = genes
    return new_genomes

def filter_otus(otus_to_keep, genomes):
    """takes a list of otus and removes all otus not in list
    """
    new_genomes = dict()
    for otu in otus_to_keep:
        if otu in genomes:
            new_genomes[otu] = genomes[otu]
    return new_genomes

def filter_by_list(genomes, list_genes):
    """reads in a list of genes separated by newlines.  filters from genomes.
    """
    f = open(list_genes, 'U')
    f = f.read()
    genes = f.strip().split('\n')
    
    return filter_genes(set(genes), genomes)
    
def top_quantile(genomes, quantile):
    """
        input: genomes = dict with genome labels as keys and a list of genes as values
               quantile = percentile of rank abundance curve below which to keep genes
        
        output: new_genomes: same as input but filtered to only contain genes that are below the
                         quantile as values and the genomes that contain them as keys
    """
    #transform to KO counter
    kos = Counter()
    for genome in genomes:
        for gene in genomes[genome]:
            kos[gene]+=1
    kos = list(kos.items())
    kos.sort(key=itemgetter(1))
    
    # #testing
    # temp = ""
    # for ko in kos:
    #     temp += str(ko[0]) + '\t' + str(ko[1]) + '\n'
    # f = open('test.txt', 'w')
    # f.write(temp)
    # f.close()
    # del temp
    # #end testing
    
    vals = [ko[1] for ko in kos]
    kos = kos[:bisect.bisect(vals, numpy.percentile(vals, quantile))]
    kos = set([ko[0] for ko in kos])
    return filter_genes(kos, genomes)

def filter_sig_OTUs(genomes, sig_file, sig_cutoff):
    """adapted from stats_parser.py from initial editing of seed code
    """
    #open file to be parsed and toss first line
    f = open(sig_file, 'U')
    if f.closed:
        "KO's not opened"
        return genomes
    f.readline()

    #for each line in file write KO, direction, weight
    otus = set()
    for line in f:
        line = line.split('\t')
        if float(line[1]) < sig_cutoff:
            otus.add(line[0])
    return filter_otus(otus, genomes)
    
def make_network(genomes, prefix):
    gg = parse_gg()
    edge_info = dict()
    nodes = dict()
    for genome, genes in genomes.iteritems():
        nodes[genome] = Node(list(genes), {'eng_name':gg[genome], 'type':'genome'})
        for gene in genes:
            edge_info[(genome, gene)] = {'type':'has'}
            if gene not in nodes:
                nodes[gene] = Node(info={'type':'gene', 'eng_name':kegg_parser.get_ko_name(gene)})
    return Network(prefix, nodes, edge_info)

def main(input_file, output_prefix, quantile, list_genes, pathway, sig_cutoff, sig_file, filt_type):
    genomes = parse_meta_contribs(input_file)
    
    print "Number of OTUs before filtering: " + str(len(genomes.keys()))
    kos = set()
    for genome in genomes:
        kos = kos | genomes[genome]
    print "Number of KOs before filtering: " + str(len(kos))
    
    if quantile != None:
        print "filtering out bottom quantile"
        genomes = top_quantile(genomes, quantile)
    
    if pathway != None:
        print "filtering to only include designated pathway"
        genomes = filter_genes(set(kegg_parser.get_kos_from_pathway(pathway[-5:])), genomes)
    
    if list_genes != None:
        print "filtering to only include specified genes"
        genomes = filter_by_list(genomes, list_genes)
    
    if filt_type == "OTU" and sig_cutoff != None and sig_file != None:
        print "filtering to only include significantly changed OTUs"
        genomes = filter_sig_OTUs(genomes, sig_file, sig_cutoff)
        
    if filt_type == "gene" and sig_cutoff != None and sig_file != None:
        print "filtering to only include significantly changed genes"
        
    
    print "Number of OTUs: " + str(len(genomes.keys()))
    kos = set()
    for genome in genomes:
        kos = genomes[genome] | kos
    print "Number of KO's: " + str(len(kos))
    make_network(genomes, output_prefix).print_network()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", help="location of metagenome contributions file")
    parser.add_argument("-o", "--output_prefix", help="prefix for output files")
    parser.add_argument("-q", "--quantile", type = float, help="quantile of rank abundance below which to keep genes")
    parser.add_argument("-g", "--list_genes", help="file with newline separated list of genes to keep")
    parser.add_argument("-O", "--list_OTUs", help="file with newline separated list of OTUs to keep")
    parser.add_argument("-p", "--pathway", help="KEGG pathway ID to create network from")
    parser.add_argument("-c", "--sig_cutoff", type = float, help="significance cutoff to include a gene/OTU")
    parser.add_argument("-f", "--sig_file", help="file containing significance values for genes/OTUs")
    parser.add_argument("-t", "--type", help="OTU or gene, to filter")
    args = parser.parse_args()
        
    main(args.input_file, args.output_prefix, args.quantile, args.list_genes, args.pathway, args.sig_cutoff, args.sig_file, args.type)