"""make_map_pathway.py"""

from biom import load_table
from make_gg_genomes import get_genome

def get_otus_from_table(biom_file):
    table = load_table(biom_file)
    return table.ids(axis = "observation")

def main(table, output):
    otus = get_otus_from_table(table)
    f = open(output, 'w')
    gene_number = 1
    for otu in otus:
        f.write("# " + otu + '\n')
        genome = get_genome(otu)
        for gene in genome:
            f.write("gene" + str(gene_number) + '\t' + gene + '\n')
            gene_number+=1
    f.close()

def one_org(table, output):
    otus = get_otus_from_table(table)
    gene_number = 1
    genes = set()
    for otu in otus:
        genome = get_genome(otu)
        genes = genes | set(genome)
    f = open(output, 'w')
    for i, gene in enumerate(genes):
        f.write("gene" + str(i) + '\t' + gene + '\n')
    f.close()