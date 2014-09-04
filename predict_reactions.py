"""predict_reactions.py
Takes in a gene biom table containing KO's and translates it into a table of KEGG reactions
"""

from biom import load_table
from biom.table import Table
from collections import defaultdict
from functools import partial
import dict_manips
import numpy
import argparse

def get_reactions():
    """get compounds for each reaction and each KO"""
    f = open("reaction", 'U')
    f = f.read()
    f = f.strip().split('///')
    pathway2ko = dict()
    ko2rxn = dict()
    rxn2co = dict()
    for entry in f:
        i = 0
        entry = entry.strip().split('\n')
        kos = []
        paths = []
        hasOrtho = False
        
        while i < len(entry):
            rev = False
            new_start = entry[i][:12].strip()
            line = entry[i][12:].strip()

            if new_start == "ENTRY":
                r = line.strip().split()[0]
                start = "ENTRY"
            elif new_start == "EQUATION":
                equ = line.split('=>')
                if equ[0][-1] == '<':
                    rev = True
                    equ[0] = equ[0][:-1]
                reacts = list()
                for part in equ[0].strip().split():
                    if part[0] == 'C' or part[0] == 'G':
                        reacts.append(part[:6])
                prods = list()
                for part in equ[1].strip().split():
                    if part[0] == 'C' or part[0] == 'G':
                        prods.append(part[:6])
                start = "EQUATION"
            elif new_start == "PATHWAY":
                paths.append(line.strip().split()[0][-5:])
                start = "PATHWAY"
            elif new_start == "ORTHOLOGY":
                kos.append(line.strip().split()[0])
                hasOrtho = True
                start = "ORTHOLOGY"
            elif new_start == "":
                if start == "ORTHOLOGY":
                    kos.append(line.strip().split()[0])
                if start == "PATHWAY":
                    paths.append(line.strip().split()[0][-5:])
            else:
                start = new_start
            i+=1
        if hasOrtho == True:
            if len(r) != 6 and r.startswith('R') == False:
                print r
            for ko in kos:
                if ko in ko2rxn:
                    ko2rxn[ko].add(r)
                else:
                    ko2rxn[ko] = set([r])
            for path in paths:
                if path in pathway2ko:
                    pathway2ko[path] = pathway2ko[path] | set(kos)
                else:
                    pathway2ko[path] = set(kos)
            rxn2co[r] = reacts,prods,rev

    return pathway2ko,ko2rxn,rxn2co
    
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

def main(input_file, output_file):
    #get reactions
    pathway2ko, ko2rxn, rxn2co = get_reactions()
    
    #parse biom table
    table = load_table(input_file)

    #combine rows of KO's from the same reaction
    table = collapse_to_reactions(table, ko2rxn)
        
    #print biom table
    table.to_json("predict_reactions.py", open(output_file, 'w'))
    
    #print to tab delimited biom table
    #f = open(output_file, 'w')
    #f.write(table.to_tsv())
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="location of gene biom table")
    parser.add_argument("-o", "--output", help="name of new reaction biom table")
    args = parser.parse_args()
    main(args.input, args.output)