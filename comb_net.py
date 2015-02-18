from make_network import metab_from_rxns_by_pathway, remove_common_cos
from biom import load_table
from network import Network

RXN_BIOM_LOC = "/Users/shafferm/lab2/asthma/biom/rarefied/metagenome/rxn_not_merged.biom"

has_rxns = load_table(RXN_BIOM_LOC).ids(axis='observation')
pah = metab_from_rxns_by_pathway('00624', has_rxns, 'pah')
nap = metab_from_rxns_by_pathway('00626', has_rxns, 'nap')
benzo = metab_from_rxns_by_pathway('00362', has_rxns, 'benzo')

comb_net = Network('comb')
comb_net.merge_network(pah)
comb_net.merge_network(nap)
comb_net.merge_network(benzo)

remove_common_cos(comb_net).print_network('output/')
