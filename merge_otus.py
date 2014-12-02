"""merge_otus.py
Take a biom table and filter to a grouping based on a list of OTU's. So for
example give a list of OTU's associated with a pathway or the metabolism of a
compound and get as output a table with the group given.
"""

from biom import load_table
from biom.table import Table
import numpy

def merged_sum(table, otus):
    table1 = table.filter(otus, axis='observation', inplace=False)
    sums = table1.sum(axis='sample')
    return sums

def main(table_loc, otu_list, collapsed_name, output_file, classic=False):
    table = load_table(table_loc)
    f = open(otu_list)
    otus = f.read().strip().split()
    otus = set(otus) & set(table.ids(axis="observation"))
    table1 = table.filter(otus, axis="observation", inplace=False)
    table2 = table.filter(otus, axis="observation", invert=True, inplace=False)
    sums1 = table1.sum(axis='sample')
    sums2 = table2.sum(axis='sample')
    new_table = Table(numpy.array([sums1,sums2]), [collapsed_name, "not_"+collapsed_name], table.ids(axis="sample"), type="otu baptable")
    
    if classic:
        # print to tab delimited biom table
        open(output_file, 'w').write(new_table.to_tsv())
    else:
        # print biom table
        new_table.to_json("predict_reactions.py", open(output_file, 'w'))
