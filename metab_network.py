"""metab_network.py
Takes group_significance.py output file with reactions and uses it to create a metabolic network 
with compounds as nodes and reactions as edges.
"""

import argparse
import parse_KEGG

def parse_group_sig(contribs_loc, otu_cat, raw_p=False):
    """
        input: sig_loc = location of group significance file
        outputs:    rxns = 
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
        rxns.add(line[0])
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
    return rxns, sig_dict

def get_list(list_):
    """get a list of genes from the file
    """
    f = open(list_)
    items = set()
    for line in f:
        items.add(line.strip())
    return items
    
def sig_changed(sigs, cutoff):
    rxns = set()
    for rxn in sigs:
        if sigs[rxn] < cutoff:
            rxns.add(rxn)
    return rxns
    
def make_network(rxns, sigs, prefix):
    """Generate a network which has compounds as nodes and edges as KO's
    """
    #get edges and compounds
    rxn2co = parse_KEGG.get_reactions()
    #rxn2co = parse_KEGG.parse_reaction_mapformula()
    rxn_names = parse_KEGG.get_rxn_names()
    #rare_cos = get_list("rare_cos2.txt")
    #rare_cos = parse_KEGG.parse_reaction_mapformula()
    edges = set()
    cos = set()
    for rxn in rxns:
        if rxn in rxn2co:
            for react in rxn2co[rxn][0]:
                #if react in rare_cos:
                    cos.add(react)
                    for prod in rxn2co[rxn][1]:
                        #if prod in rare_cos:
                            cos.add(prod)
                            if rxn2co[rxn][2] == True:
                                edges.add((react, prod, rxn, rxn_names[rxn], str(sigs[rxn]), "true"))
                            else:
                                edges.add((react, prod, rxn, rxn_names[rxn], str(sigs[rxn]), "false"))
    
    #create and print edges file
    new_edges = list()
    for edge in edges:
        new_edges.append('\t'.join(edge))
    
    print "Number of edges: " + str(len(new_edges))
    
    f = open(prefix+"_edges.txt", 'w')
    f.write("reactant\tproduct\trxn_id\treaction_name\tFDR\treversible\n" + '\n'.join(new_edges) + '\n')
    f.close()
    
    #create and print nodes file
    nodes = set()
    co_info = parse_KEGG.get_co_info()
    for co in cos:
        try:
            nodes.add((co, co_info[co].name))
        except KeyError:
            nodes.add((co, co))
    
    new_nodes = list()
    for node in nodes:
        new_nodes.append('\t'.join(node))
    
    print "Number of nodes: " + str(len(new_nodes))
        
    f = open(prefix+"_nodes.txt", 'w')
    f.write("compound\tcompound_name\n" + '\n'.join(new_nodes) + '\n')
    f.close()

def main(input_file, output_prefix, rxn_list, pathway, sig_cutoff, otu_cat, raw_p):
    #setup: parse reactions, parse out kos from meta_contribs file
    rxns, sigs = parse_group_sig(input_file, otu_cat, raw_p)
    
    print "Number of rxns before filtering: " + str(len(rxns))
    
    #filter by pathway
    if pathway != None:
        rxns = rxns & parse_KEGG.get_pathway2rxns()[pathway]
        
    #filter by list
    if rxn_list != None:
        rxns = rxns & get_list(rxn_list)
    
    if sig_cutoff != None:
        rxns = rxns & sig_changed(sigs, sig_cutoff)
        
    print "Number of reactions: " + str(len(rxns))
    make_network(rxns, sigs, output_prefix)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="location of group_significance.py output containing reactions", required=True)
    parser.add_argument("-o", "--output_prefix", help="prefix for output files", required=True)
    parser.add_argument("-l", "--list_rxns", help="file with newline separated list of rxns to keep")
    parser.add_argument("-p", "--pathway", help="KEGG pathway ID to create network from")
    parser.add_argument("-c", "--sig_cutoff", type = float, help="significance cutoff to include a reaction")
    parser.add_argument("--otu_cat", action='store_true', help="if data came from otu_category_significance.py output files")
    parser.add_argument("--raw_p", action='store_true', help="if the raw p-value should be used instead of the FDR corrected value")
    args = parser.parse_args()
    main(args.input, args.output_prefix, args.list_rxns, args.pathway, args.sig_cutoff, args.otu_cat, args.raw_p)