"""anal_comp_metab.py"""

from parse_KEGG import KEGG_Parser
from make_gg_genomes import get_genome
from biom import load_table

kegg_parser = KEGG_Parser()

def get_genomes(otu_list):
    genomes = dict()
    for otu in otu_list:
        genomes[otu] = get_genome(otu)
    return genomes

def get_first_rxn(genome, compound):
    comp_genes = set()
    for rxn in kegg_parser.get_co_info(compound).rxn:
        comp_genes = comp_genes | kegg_parser.get_kos_from_rxn(rxn)
    first_genes = set()
    for gene in genome:
        if gene in comp_genes:
            first_genes.add(gene)
    return first_genes

def can_react_with(genomes, compound):
    can_otus = list()
    can_kos = dict()
    for genome in genomes:
        if len(get_first_rxn(genomes[genome], compound))>0:
            can_kos[genome] = len(get_first_rxn(genomes[genome], compound))
            can_otus.append(genome)
    return filter_by_list(genomes, can_otus), can_kos
    
def has_pathway_percent(genomes, pathway, cutoff):
    pathway = kegg_parser.get_kos_from_pathway(pathway)
    otus_to_keep = list()
    for genome in genomes:
        percent = float(len(set(genomes[genome])&pathway))/float(len(pathway))
        if percent > cutoff:
            otus_to_keep.append(genome)
            print str(percent) + '\t' + genome
    return filter_by_list(genomes, otus_to_keep)

def filter_by_list(genomes, otus_to_keep):
    return {otu:genes for (otu,genes) in genomes.iteritems() if otu in otus_to_keep}

def get_otus_from_table(biom_file):
    table = load_table(biom_file)
    return table.ids(axis = "observation")

def filter_kos(genomes, kos_file):
    kos_to_toss = set(open(kos_file, 'U').read().split())
    for genome in genomes:
        genomes[genome] = list(set(genomes[genome])-kos_to_toss)
    return genomes

def main(biom_file, output, compound, pathway, cutoff=0., kos_to_filter=None):
    # get otus from biom table
    otu_list = get_otus_from_table(biom_file)
    print "retrieved otus"
    # get genomes for otus
    genomes = get_genomes(otu_list)
    print "retrieved genomes for otus: " + str(len(genomes))
    # filter out undesired kos
    if kos_to_filter != None:
        genomes = filter_kos(genomes, kos_to_filter)
        print "filtered out KO's"
    # filter out genomes which can't react with compound
    genomes = can_react_with(genomes, compound)
    print "filtered to genomes that can react with " + compound + ": " + str(len(genomes))
    # filter out genomes with pathways with under a given percent of genes
    genomes = has_pathway_percent(genomes, pathway, cutoff)
    print "filtered to genomes that contain cutoff percent of genes from pathway: " + str(len(genomes))
    open(output, 'w').write('\n'.join(genomes.keys())+'\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True,
                        help="location of input biom table")
    parser.add_argument("-o", "")
    parser.add_argument("-c", "--coumpound", required=True,
                        help="KEGG compound id")
    parser.add_argument("-p", "--pathway", required=True,
                        help="KEGG pathway number")
    parser.add_argument("--cutoff", default=.6,
                        help="percent pathway presensce required to keep")
    parser.add_argument("--filter_kos",
                        help="file containing list of KO's to filter separated by whitespace")
    
    main(args.input, args.output, args.compound, args.pathway, args.cutoff, args.filter_frequent)
