"""make_network.py"""

import parse_KEGG as k

kegg_parser = k.KEGG_Parser()

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

def metab_from_genes(all_genes, genes_to_keep, prefix):
    # node header format: source (tab) target (tab) if reaction (tab) if present (tab) name
    node_header = ("name", "reaction", "present", "comp_name")
    rxns = set()
    for gene in all_genes:
        rxns = rxns | kegg_parser.get_rxns_from_ko(gene)
    has_rxns = set()
    for gene in genes_to_keep:
        has_rxns = has_rxns | set(kegg_parser.get_rxns_from_ko(gene))
    
    edges = set()
    nodes = set()
    for rxn in rxns:
        reacts, prods, rev = kegg_parser.get_rxn(rxn)
        for react in reacts:
            edges.add((react, rxn))
            nodes.add((react, False, False, kegg_parser.get_co_info(react).name))
        for prod in prods:
            edges.add((rxn, prod))
            nodes.add((prod, False, False, kegg_parser.get_co_info(prod).name))
        if rxn in has_rxns:
            nodes.add((rxn, True, True, kegg_parser.get_rxn_name(rxn)))
        else:
            nodes.add((rxn, False, True, kegg_parser.get_rxn_name(rxn)))
    print_network(prefix, edges, nodes, node_header=node_header)

def metab_from_genes_by_pathway(pathway, genes, prefix):
    # node header format: source (tab) target (tab) if reaction (tab) if present (tab) name
    node_header = ("name", "reaction", "present", "comp_name")
    rxns = kegg_parser.get_rxns_from_pathway(pathway)
    has_rxns = set()
    for gene in genes:
        gene_rxns = kegg_parser.get_rxns_from_ko(gene)
        has_rxns = has_rxns | gene_rxns
    edges = set()
    nodes = set()
    for rxn in rxns:
        reacts, prods, rev = kegg_parser.get_rxn(rxn)
        for react in reacts:
            edges.add(list(react, rxn))
            nodes.add(list(react, False, False, kegg_parser.get_co_info(react).name))
        for prod in prods:
            edges.add(list(rxn, prod))
            nodes.add(list(prod, False, False, kegg_parser.get_co_info(prod).name))
        if rxn in has_rxns:
            nodes.add(list(rxn, True, True, kegg_parser.get_rxn_name(rxn)))
        else:
            nodes.add(list(rxn, False, True, kegg_parser.get_rxn_name(rxn)))
    print_network(prefix, edges, nodes, node_header=node_header)

def add_sigs_to_network(sig_file, nodes, node_header=None):
    sig_dict = parse_group_sig(sig_file)
    for node in nodes:
        node.append(sig_dict[node[0]])
    if node_header == None:
        return nodes
    else:
        node_header.append("sig_level")
        return nodes, node_header

def print_network(prefix, edges, nodes, edge_header=None, node_header=None):
    # edges
    f = open(prefix+'_edges.txt', 'w')
    new_edges= []
    if edge_header != None:
        new_edges.apped('\t'.join(edge_header))
    for edge in edges:
        edge = map(str, edge)
        new_edges.append('\t'.join(edge))
    f.write('\n'.join(new_edges) + '\n')
    f.close()

    # nodes
    f = open(prefix+'_nodes.txt', 'w')
    new_nodes = []
    if node_header != None:
        new_nodes.append('\t'.join(node_header))
    for node in nodes:
        node = map(str, node)
        new_nodes.append('\t'.join(node))
    f.write('\n'.join(new_nodes) + '\n')
    f.close()