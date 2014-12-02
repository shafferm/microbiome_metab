def merged_test():
    from anal_comp_metab import *
    from merge_otus import merged_sum
    from biom import load_table
    from predict_reactions import make_biom_table
    
    table = load_table("/Users/shafferm/lab2/asthma/biom/rarefied/otu_closed_rar_norm_s4_n2.biom")
    table_otus = table.ids(axis='observation')
    genomes = get_genomes(table_otus)

    comps = ['C00829', 'C14222', 'C07535', 'C14315', 'C11422']
    comp_otus = dict()
    for comp in comps:
        comp_otus[comp] = can_react_with(genomes, comp).keys()
    comp_sums = dict()
    for comp, otus in comp_otus.iteritems():
        comp_sums[comp] = merged_sum(table, otus)
        print comp_sums[comp]
    new_table = make_biom_table(comp_sums, table.ids())
    new_table.to_json("merge_sums", open("test_sums.biom", 'w'))

def net_test():
    import network as n
    a = n.Node(['B'])
    b = n.Node(['A'])
    net = n.Network('basic')
    net.add_node('A', a)
    net.add_node('B', b)
    net.print_network()
