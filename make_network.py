"""make_network.py"""

import parse_KEGG as k
from collections import defaultdict

kegg_parser = k.KEGG_Parser()
from network import Network, Node

COMMON_COS_LOC = "common_cos.txt"

def parse_group_sig(contribs_loc, otu_cat=False, raw_p=False):
    """
        input: sig_loc = location of group significance file
        outputs:    sig_dict = 
    """
    f = open(contribs_loc, 'U')
    if f.closed == True:
        print "file location does not exist"
        return set()
    f.readline() #toss header
    rxns = set()
    sig_dict = dict()
    for line in f:
        line = line.strip().split('\t')
        if otu_cat == False:
            #FDR
            if raw_p == False:
                sig_dict[line[0]] = float(line[3])
            #raw p
            else:
                sig_dict[line[0]] = float(line[2])
        else:
            if raw_p == False:
                #FDR
                sig_dict[line[0]] = float(line[6])
            else:
                #raw p
                sig_dict[line[0]] = float(line[1])
    f.close()
    return sig_dict

def metab_from_rxns(rxns, prefix='net', has_rxns=None, rxns_as_nodes=True):
    nodes = dict()
    if rxns_as_nodes == True:
        for rxn in rxns:
            reacts, prods, rev = kegg_parser.get_rxn(rxn)
            nodes[rxn] = Node()
            nodes[rxn].add_info('type', 'reaction')
            nodes[rxn].add_info('eng_name', kegg_parser.get_rxn_name(rxn))
            if has_rxns != None:
                if rxn in has_rxns:
                    nodes[rxn].add_info('present', True)
                else:
                    nodes[rxn].add_info('present', False)
            for react in reacts:
                if react in nodes:
                    nodes[react].add_target(rxn)
                else:
                    nodes[react] = Node([rxn])
                nodes[react].add_info('type', 'compound')
                nodes[react].add_info('present', False)
                nodes[react].add_info('eng_name', kegg_parser.get_co_info(react).name)
            for prod in prods:
                nodes[rxn].add_target(prod)
                if prod not in nodes:
                    nodes[prod] = Node()
                nodes[prod].add_info('type', 'compound')
                nodes[prod].add_info('present', False)
                nodes[prod].add_info('eng_name', kegg_parser.get_co_info(prod).name)
    else:
        edge_info = defaultdict(dict)
        for rxn in rxns:
            reacts, prods, rev = kegg_parser.get_rxn(rxn)
            for react in reacts:
                if react not in nodes:
                    nodes[react] = Node()
                    nodes[react].add_info('eng_name', kegg_parser.get_co_info(react).name)
                for prod in prods:
                    nodes[react].add_target(prod)
                    edge_info[(react, prod)]['rxn_id'] = rxn
                    edge_info[(react, prod)]['eng_name'] = kegg_parser.get_rxn_name(rxn)
                    if has_rxns != None:
                        edge_info[(react, prod)]['present'] = rxn in has_rxns
                    if prod not in nodes:
                        nodes[prod] = Node()
                    nodes[prod].add_info('eng_name', kegg_parser.get_co_info(prod).name)
    return Network(prefix, nodes, edge_info)

def metab_from_rxns_co_only(rxns, has_rxns=None, prefix='net'):
    nodes = dict()
    edge_info = dict()
    for rxn in rxns:
        reacts, prods, rev = kegg_parser.get_rxn(rxn)
        for react in reacts:
            if react not in nodes:
                nodes[react] = Node()
                nodes[react].add_info('eng_name', kegg_parser.get_co_info(react).name)
            for prod in prods:
                nodes[react].add_target(prod)
                edge_info[(react, prod)]['rxn_id'] = rxn
                edge_info[(react, prod)]['eng_name'] = kegg_parser.get_rxn_name(rxn)
                edge_info[(react, prod)]['present'] = rxn in has_rxns
                if prod not in nodes:
                    nodes[prod] = Node()
                nodes[prod].add_info('eng_name', kegg_parser.get_co_info(prod).name)
    return Network(prefix, nodes)

def metab_from_rxns_by_pathway(pathway, has_rxns=None, prefix='net', rxns_as_nodes=True):
    rxns = kegg_parser.get_rxns_from_pathway(pathway)
    return metab_from_rxns(rxns, has_rxns=has_rxns, prefix=prefix, rxns_as_nodes=rxns_as_nodes)

def metab_from_genes(all_genes, genes_to_keep, prefix):
    # node header format: source (tab) target (tab) if reaction (tab) if present (tab) name
    rxns = set()
    for gene in all_genes:
        rxns = rxns | kegg_parser.get_rxns_from_ko(gene)
    has_rxns = set()
    for gene in genes_to_keep:
        has_rxns = has_rxns | set(kegg_parser.get_rxns_from_ko(gene))
    return metab_from_rxns(rxns, has_rxns=has_rxns, prefix=prefix)

def metab_from_genes_by_pathway(pathway, genes, prefix):
    has_rxns = set()
    for gene in genes:
        gene_rxns = kegg_parser.get_rxns_from_ko(gene)
        has_rxns = has_rxns | gene_rxns
    return metab_from_rxns_by_pathway(pathway, has_rxns=has_rxns, prefix=prefix)

def add_sigs_to_nodes(sig_file, network, otu_cat=False, raw_p=False):
    sig_dict = parse_group_sig(sig_file, otu_cat, raw_p)
    for name, node in network.nodes.iteritems():
        try:
            node.add_info('sig', sig_dict[name])
        except:
            pass
    return network

def remove_common_cos(network):
    common_cos = open(COMMON_COS_LOC).read().strip().split()
    network.remove_nodes(common_cos)
    return network
