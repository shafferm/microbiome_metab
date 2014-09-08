"""metab_network.py
Takes metagenome_contributions.py output file from PICRUSt and uses it to create a metabolic
network with compounds as nodes and reactions as edges.
"""

import argparse
import parse_reaction

def parse_meta_contribs(contribs_loc):
    """
        input: contribs_loc = location of metagenome contributions file
        outputs:    genomes = dict with OTU's as keys and sets of ko's as values
    """
    f = open(contribs_loc, 'U')
    if f.closed == True:
        print "file location does not exist"
        return set()
    f.readline() #toss header
    kos = set()
    for line in f:
        line = line.strip().split('\t')
        kos.add(line[0])
    f.close()
    return kos

def filter_genes(genes_to_keep, genomes):
    """takes a list of genes to keep and removes all others from the genome, takes in a set of
    genes so the input gene data structure must be cast to a set.
    """
    new_genes = Counter()
    for genome in genomes:
        genes = genes_to_keep & genomes[genome]
        if len(genes) > 0:
            new_genomes[genome] = genes
    return new_genomes

def list_genes(gene_list):
    """get a list of genes from the file
    """
    f = open(gene_list)
    genes = set()
    for line in f:
        genes.add(line.strip())
    return genes
    
def sig_changed(sig_file, cutoff):
    kos = set()
    f = open(sig_file, 'U')
    f.readline()

    #for each line in file write KO, direction, weight
    for line in f:
        line = line.split('\t')
        if float(line[1]) < sig_cutoff:
            kos.add(line[0])
    return kos
    
def make_network(kos):
    """Generate a network which has compounds as nodes and edges as KO's
    """
    ko2rxn = parse_reaction.get_ko2rxns()
    rxn2co = parse_reaction.get_reactions()
    edges = set()
    for ko in kos:
        if ko in ko2rxn:
            for rxn in ko2rxn[ko]:
                if rxn in rxn2co:
                    for react in rxn2co[rxn][0]:
                        for prod in rxn2co[rxn][1]:
                            edges.add((react, prod))
                            if rxn2co[rxn][2] == True:
                                edges.add((prod, react))
    new_edges = list()
    for edge in edges:
        new_edges.append('\t'.join(edge))
        
    new_edges = '\n'.join(new_edges)
    f = open("edges.txt", 'w')
    f.write(new_edges+'\n')
    f.close()

def main(input_file, output_prefix, gene_list, pathway, sig_cutoff):
    #setup: parse reactions, parse out kos from meta_contribs file
    kos = parse_meta_contribs(input_file)
    
    print "Number of KOs before filtering: " + str(len(kos))
    
    #filter by pathway
    if pathway != None:
        kos = kos & parse_reaction.get_pathway2kos()[pathway]
        
    #filter by list
    if gene_list != None:
        kos = kos & list_genes(gene_list)
    
    if sig_cutoff != None:
        kos = kos & sig_changed(input_sig, sig_cutoff)
        
    print "Number of KO's: " + str(len(kos))
    make_network(kos)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="location of metagenome contributions file or file containing significance values for genes/OTUs")
    parser.add_argument("-o", "--output_prefix", help="prefix for output files")
    parser.add_argument("-l", "--list_genes", help="file with newline separated list of genes to keep")
    parser.add_argument("-p", "--pathway", help="KEGG pathway ID to create network from")
    parser.add_argument("-c", "--sig_cutoff", type = float, help="significance cutoff to include a gene")
    args = parser.parse_args()
    main(args.input, args.output_prefix, args.list_genes, args.pathway, args.sig_cutoff)