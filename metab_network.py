"""metab_network.py
Takes group_significance.py output file with reactions and uses it to create a metabolic network 
with compounds as nodes and reactions as edges.
"""

import argparse
import parse_reaction

def parse_group_sig(contribs_loc):
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
        sig_dict[line[0]] = float(line[3])
    f.close()
    return rxns, sig_dict

def list_rxns(rxn_list):
    """get a list of genes from the file
    """
    f = open(gene_list)
    rxns = set()
    for line in f:
        rxns.add(line.strip())
    return rxns
    
def sig_changed(sigs, cutoff):
    rxns = set()
    for rxn in sigs:
        if sigs[rxn] < cutoff:
            rxns.add(rxn)
    return rxns
    
def make_network(rxns, sigs, prefix):
    """Generate a network which has compounds as nodes and edges as KO's
    """
    rxn2co = parse_reaction.get_reactions()
    edges = set()
    for rxn in rxns:
        if rxn in rxn2co:
            for react in rxn2co[rxn][0]:
                for prod in rxn2co[rxn][1]:
                    if rxn2co[rxn][2] == True:
                        edges.add((react, prod, rxn, str(sigs[rxn]), "true"))
                    else:
                        edges.add((react, prod, rxn, str(sigs[rxn]), "false"))
    
    new_edges = list()
    for edge in edges:
        new_edges.append('\t'.join(edge))
        
    f = open(prefix+"_edges.txt", 'w')
    f.write("reactant\tproduct\trxn_id\tFDR\treversible\n" + '\n'.join(new_edges) + '\n')
    f.close()

def main(input_file, output_prefix, rxn_list, pathway, sig_cutoff):
    #setup: parse reactions, parse out kos from meta_contribs file
    rxns, sigs = parse_group_sig(input_file)
    
    print "Number of rxns before filtering: " + str(len(rxns))
    
    #filter by pathway
    if pathway != None:
        kos = kos & parse_reaction.get_pathway2rxns()[pathway]
        
    #filter by list
    if rxn_list != None:
        rxns = rxns & list_genes(gene_list)
    
    if sig_cutoff != None:
        rxns = rxns & sig_changed(sigs, sig_cutoff)
        
    print "Number of reactions: " + str(len(rxns))
    make_network(rxns, sigs, output_prefix)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="location of group_significance.py output containing reactions")
    parser.add_argument("-o", "--output_prefix", help="prefix for output files")
    parser.add_argument("-l", "--list_rxns", help="file with newline separated list of rxns to keep")
    parser.add_argument("-p", "--pathway", help="KEGG pathway ID to create network from")
    parser.add_argument("-c", "--sig_cutoff", type = float, help="significance cutoff to include a reaction")
    args = parser.parse_args()
    main(args.input, args.output_prefix, args.list_rxns, args.pathway, args.sig_cutoff)