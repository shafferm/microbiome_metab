import otu_ko_network as ok
import make_network as mn
from biom import load_table
from anal_comp_metab import get_genomes

OTU_BIOM = "/Users/shafferm/lab2/asthma/biom/rarefied/otu_closed_rar_norm.biom"
RXN_BIOM = "/Users/shafferm/lab2/asthma/biom/rarefied/metagenome/rxn_s4_n2.biom"

genome_list = load_table(OTU_BIOM).ids(axis="observation")
genomes = get_genomes(genome_list)
otu_ko_network = ok.make_network(genome, "otu-ko")

rxn_list = load_table(RXN_BIOM).ids(axis="observation")
rxn_network = mn.metab_from_rxns(rxn_list, prefix="metab")


