import otu_ko_network as ok
import make_network as mn
from biom import load_table
from anal_comp_metab import get_genomes
from parse_KEGG import KEGG_Parser
from network import Network

OTU_BIOM = "/Users/shafferm/lab2/asthma/biom/rarefied/otu_closed_rar_norm.biom"
RXN_BIOM = "/Users/shafferm/lab2/asthma/biom/rarefied/metagenome/rxn_not_merged.biom"

kegg_parser = KEGG_Parser()

genome_list = load_table(OTU_BIOM).ids(axis="observation")
genomes = get_genomes(genome_list)
otu_ko_net = ok.make_network(genomes, "otu-ko")



rxn_list = load_table(RXN_BIOM).ids(axis="observation")
rxn_net = mn.metab_from_rxns(rxn_list, prefix="metab")

mega_net = Network('mega_net', nodes=otu_ko_net.nodes, edge_info=otu_ko_net.edge_info)
mega_net.add_node_dict(rxn_net.nodes)

for name, node in mega_net.nodes.iteritems():
    if node.info['type'] == 'gene':
        node_rxns = kegg_parser.get_rxns_from_ko(name)
        for rxn in node_rxns:
            if rxn in mega_net.nodes:
                mega_net.nodes[name].add_target(rxn)

mega_net.print_network()
filt_mega_net.print_network()

#otu_ko.print_network()
#filt_otu_ko.print_network()

#metab.print_network()
#filt_metab.print_network()