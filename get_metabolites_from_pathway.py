"""get_metabolites_from_pathway.py
"""

from parse_KEGG import get_pathway2kos, get_ko2rxns, get_reactions, get_co_info
import argparse

def main(pathways, output):
    pathway2kos = get_pathway2kos()
    kos = set()
    for pathway in pathways:
        kos = kos | set(pathway2kos[pathway])
    del pathway2kos
    ko2rxns = get_ko2rxns()
    rxns = set()
    for ko in kos:
        rxns = rxns | set(ko2rxns[ko])
    del ko2rxns
    reactions = get_reactions()
    metabolites = set()
    for rxn in rxns:
        reacts, prods, rever = reactions[rxn]
        metabolites = metabolites | set(reacts) | set(prods)
    del reactions
    metab_info = ['KEGG_id\tname\tformula\tmass']
    co_info = get_co_info()
    for metabolite in metabolites:
        compound = co_info[metabolite]
        metab_info.append('\t'.join([str(compound.co), str(compound.name), str(compound.formula), str(compound.mass)]))
    f = open(output, 'w')
    f.write('\n'.join(metab_info)+'\n')
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--pathways", nargs='+', required=True,
                        help="pathways to get metabolites from")
    parser.add_argument("-o", "--output", required=True,
                        help="file to output results")
    args = parser.parse_args()
    
    main(args.pathways, args.output)