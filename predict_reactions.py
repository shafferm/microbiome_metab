"""predict_reactions.py
Takes in a gene biom table containing KO's and translates it into a table of KEGG reactions

TODO: refactor to use get_reactions.py from parse_reaction.py
"""

from biom import load_table
from biom.table import Table
from collections import defaultdict
from functools import partial
import dict_manips
import numpy
import argparse
import parse_reaction
    
def collapse_to_reactions(table, ko2rxns):
    #get dictionary of rxns to kos
    rxns = defaultdict(partial(numpy.array, numpy.zeros(len(table.ids()))))
    for values, ko, metadata in table.iter(axis='observation'):
        if str(ko) in ko2rxns:
            for rxn in ko2rxns[str(ko)]:
                rxns[rxn] += values
    
    #merge reactions which have the same ko set
    kos2rxns = defaultdict(list)
    rxn2ko = dict_manips.reverse_dict_of_iterables(ko2rxns)
    for rxn in rxns:
        kos2rxns[frozenset(rxn2ko[rxn])].append(rxn)
    rxn_sets = kos2rxns.values()
    merged_rxns = dict()
    for rxn in rxn_sets:
        if numpy.sum(rxns[rxn[0]]) > 0:
            merged_rxns[','.join(rxn)] = rxns[rxn[0]]
    print len(merged_rxns)
    
    #create biom table
    merged_rxns = merged_rxns.items()
    data = numpy.array([x[1] for x in merged_rxns], dtype = float)
    return Table(data, [x[0] for x in merged_rxns], table.ids(), type="OTU table")

def main(input_file, output_file, classic):
    #parse biom table
    table = load_table(input_file)

    #combine rows of KO's from the same reaction
    table = collapse_to_reactions(table, parse_reaction.get_ko2rxns())
    
    if classic:
        #print to tab delimited biom table
        f = open(output_file, 'w')
        f.write(table.to_tsv())
    else:
        #print biom table
        table.to_json("predict_reactions.py", open(output_file, 'w'))
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="location of gene biom table")
    parser.add_argument("-o", "--output", help="name of new reaction biom table")
    parser.add_argument("-b", "--to_classic_table", action="store_true", 
        help="output as a tab delimited table")
    args = parser.parse_args()
    main(args.input, args.output, args.to_classic_table)