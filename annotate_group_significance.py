"""annotate_group_significance.py
Adds information to group_significance.py or otu_category_significance.py script output.
"""

import parse_reaction
import sys
import argparse

def parse_group_significance(file_loc):
    """"""
    f = open(file_loc, 'U')
    header = f.readline().strip().split('\t')
    entries = list()
    for line in f:
        line = line.strip().split('\t')
        line[0] = line[0].split(',')
        entries.append(line)
    return header, entries
    
def annotate_rxns(header, entries):
    """"""
    rxn_names = parse_reaction.get_rxn_names()
    header.append("reaction_name")
    new_entries = list()
    for entry in entries:
        names = list()
        for rxn in entry[0]:
            names.append(rxn_names[rxn])
        names = ", ".join(names)
        entry[0] = ','.join(entry[0])
        entry.append(names)
        new_entries.append(entry)
    return header, new_entries
    
def annotate_kos(header, entries):
    """"""
    ko_names = parse_reaction.get_ko_names()
    header.append("ko_name")
    new_entries = list()
    for entry in entries:
        names = list()
        for ko in entry[0]:
            names.append(ko_names[ko])
        names = ", ".join(names)
        entry[0] = ','.join(entry[0])
        entry.append(names)
        new_entries.append(entry)
    return header, new_entries

def print_annotated_file(header, entries, output):
    """EACH ENTRY MUST CONTAIN ONLY STRINGS"""
    lines = list()
    for entry in entries:
        lines.append('\t'.join(entry))
    header = '\t'.join(header) + '\n'
    lines = '\n'.join(lines) + '\n'
    lines = header + lines
    f = open(output, 'w')
    f.write(lines)

def main(in_file, out_file, kind):
    #do things
    header, entries = parse_group_significance(in_file)
    
    if kind == "reaction":
        header, entries = annotate_rxns(header, entries)
    elif kind == "compound":
        header, entries = annotate_compounds(header, entries)
    elif kind == "KO":
        header, entries = annotate_KOs(header, entries)
    else:
        print kind + " is not a valid type"
        sys.exit()
        
    print_annotated_file(header, entries, out_file)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("kind", help="type of things being compared")
    parser.add_argument("-i", "--input", help="location of group_significance.py results")
    parser.add_argument("-o", "--output", help="name of new annoated results")
    args = parser.parse_args()
    main(args.input, args.output, args.kind)