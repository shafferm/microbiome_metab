"""filter_KOs_by_pathway.py
takes a biom table and filteres it based on the KO's present
"""

from biom import load_table
from biom import Table
from biom.util import biom_open
from parse_KEGG import get_pathway2kos
import argparse
import numpy as np


def main(table_in, table_out, pathways, to_classic):
    # setup
    table = load_table(table_in)
    pathway_dict = get_pathway2kos()

    # get set of kos from pathways
    pathways_kos = set()
    for pathway in pathways:
        pathways_kos = pathways_kos | pathway_dict[pathway.strip()[-5:]]

    # get selected kos
    kos_to_keep = set(table.ids('observation')) & \
        pathways_kos
    if len(kos_to_keep) == 0:
        raise EmptySetERROR('Intersection created empty set')
    obs_ids = np.array(list(kos_to_keep))
    data = np.empty([len(obs_ids), len(table.ids('sample'))])
    for i, obs in enumerate(obs_ids):
        data[i] = table.data(obs, 'observation')

    # output
    new_table = Table(data, obs_ids, table.ids('sample'), type="OTU table")
    if to_classic:
        # print to tab delimited biom table
        f = open(table_out, 'w')
        f.write(new_table.to_tsv())
    else:
        # print json biom table
        new_table.to_json("filter_KOs_by_pathway.py", open(table_out, 'w'))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True,
                        help="location of gene biom table")
    parser.add_argument("-o", "--output", required=True,
                        help="name of new reaction biom table")
    parser.add_argument("-p", "--pathways", nargs='+', required=True,
                        help="pathway id of pathway to keep")
    parser.add_argument("-b", "--to_classic_table", action="store_true",
                        help="output as a tab delimited table")
    args = parser.parse_args()
    main(args.input, args.output, args.pathways, args.to_classic_table)
