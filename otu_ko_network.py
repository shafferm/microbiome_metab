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

def get_reactions():
    """get compounds for each reaction and each KO"""
    f = open("reaction", 'U')
    f = f.read()
    f = f.strip().split('///')
    pathway2ko = dict()
    ko2rxn = dict()
    rxn2co = dict()
    for entry in f:
        i = 0
        entry = entry.strip().split('\n')
        kos = []
        paths = []
        hasOrtho = False
        
        while i < len(entry):
            rev = False
            new_start = entry[i][:12].strip()
            line = entry[i][12:].strip()

            if new_start == "ENTRY":
                r = line.strip().split()[0]
                start = "ENTRY"
            elif new_start == "EQUATION":
                equ = line.split('=>')
                if equ[0][-1] == '<':
                    rev = True
                    equ[0] = equ[0][:-1]
                reacts = list()
                for part in equ[0].strip().split():
                    if part[0] == 'C' or part[0] == 'G':
                        reacts.append(part[:6])
                prods = list()
                for part in equ[1].strip().split():
                    if part[0] == 'C' or part[0] == 'G':
                        prods.append(part[:6])
                start = "EQUATION"
            elif new_start == "PATHWAY":
                paths.append(line.strip().split()[0][-5:])
                start = "PATHWAY"
            elif new_start == "ORTHOLOGY":
                kos.append(line.strip().split()[0][-5:])
                hasOrtho = True
                start = "ORTHOLOGY"
            elif new_start == "":
                if start == "ORTHOLOGY":
                    kos.append(line.strip().split()[0][-5:])
                if start == "PATHWAY":
                    paths.append(line.strip().split()[0][-5:])
            else:
                start = new_start
            i+=1
        if hasOrtho == True:
            for ko in kos:
                if ko in ko2rxn:
                    ko2rxn[ko].add(r)
                else:
                    ko2rxn[ko] = set(r)
            for path in paths:
                if path in pathway2ko:
                    pathway2ko[path] = pathway2ko[path] | set(kos)
                else:
                    pathway2ko[path] = set(kos)
            rxn2co[r] = reacts,prods,rev

    return pathway2ko,ko2rxn,rxn2co
    
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
    del vals
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
    
def make_network(genomes, prefix, ko2rxn, rxn2co):
    """
        input: genomes = dict with genomes as keys and sets of KO's as values 
               prefix = prefix for output files
        output: None
    """
    edges = open(prefix+"_edges.txt", 'w')
    edges.write("node1\ttype\tnode2\n")
    nodes = open(prefix+"_nodes.txt", 'w')
    nodes.write("id\ttype\n")
    nodes_added = set()
    rxns_added = set()
    for genome in genomes:
        nodes.write(genome + '\tOTU\n')
        for gene in genomes[genome]:
            edges.write(genome + '\thas\t' + gene + '\n')
            if gene not in nodes_added:
                nodes.write(gene + '\tKO\n')
                nodes_added.add(gene)
                for rxn in ko2rxn[gene]:
                    if rxn in rxn2co and rxn not in rxns_added:
                        for co in rxn2co[rxn][0]:
                            edges.write(co + '\t>\t' + gene + '\n') #reactant to gene
                        for co in rxn2co[rxn][1]:
                            edges.write(gene + '\t>\t' + co + '\n') #gene to reactant
                        if rxn2co[rxn][2]:
                            for co in rxn2co[rxn][0]:
                                edges.write(gene + '\t>\t' + co + '\n') #reverse reactant to gene
                            for co in rxn2co[rxn][1]:
                                edges.write(co + '\t>\t' + gene + '\n') #reverse gene to product
                        rxns_added.add(rxn)
    edges.close()
    nodes.close()

def main(input_file, output_prefix, quantile, list_genes, pathway, sig_cutoff, sig_file, filt_type):
    pathways, ko2rxn, rxn2co = get_reactions()
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
        genomes = filter_genes(set(pathways[pathway[-5:]]), genomes)
    
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
    make_network(genomes, output_prefix, ko2rxn, rxn2co)

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
    
#deprecated
def get_ko2rxns():
    """read in pathway information from ko2rn.xl"""
    f = open("ko2rn.xl", 'U')
    ko2rxn = dict()
    f.readline()
    for line in f:
        line = line.strip().split('\t')
        ko = line[0]
        rxns = line[1][4:-1].split()
        ko2rxn[ko] = rxns
    return ko2rxn
    
def get_rxn2cos():
    """read in CO ID's associated with each reaction"""
    rxn2co = dict()
    with open("br08202.keg", 'U') as f:
        for line in f:
            if line[0] == "D":
                line = line.strip().split('\t')
                if len(line) > 2:
                    rxn = line[0].split()[1]
                    cos = line[1].split()[0], line[2].split()[0]
                    rxn2co[rxn] = cos
    return rxn2co

def get_pathways():
    """reads in pathway information containing all KO's associated with the pathway from 
    ko00001.keg.  Location is currently assumed to be in the same folder as the script.
    """
    pathways = dict()
    with open("ko00001.keg", 'U') as f:
        for line in f:
            if line[0] == "C":
                pathway = line.split()[1]
            if line[0] == "D":
                gene = line.split()[1]
                if pathway in pathways:
                    pathways[pathway].append(gene)
                else:
                    pathways[pathway] = [gene]
    return pathways