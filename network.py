"""network.py"""

class Network:
    """class to hold network information include nodes, edge info and name.
    nodes is a dict of node objects
    edge info keys are tuples used to annotate edges as defined implicitly by the nodes
    TODO: Make nodes into Panada dataframes
    """
    
    def __init__(self, name, nodes=dict(), edge_info=dict()):
        self.name = name
        self.nodes = nodes
        self.edge_info = edge_info
            
    def generate_node_header(self):
        header = set()
        for name, node in self.nodes.iteritems():
            header = header | set(node.info.keys())
        return header
    
    def generate_edge_header(self):
        """generate edge header based on edge_info dict"""
        header = set()
        for edge in self.edge_info:
            header = header | set(self.edge_info[edge].keys())
        return header
    
    def generate_edges(self):
        edges = set()
        for name, node in self.nodes.iteritems():
            for target in node.targets:
                edges.add((name, target))
        return edges
    
    def add_node_dict(self, node_dict):
        for key, value in node_dict.iteritems():
            self.add_node(key, value)
    
    def add_node(self, name, info):
        self.nodes[name] = info
    
    def print_network(self):
        edge_names = set()
        nodes = list()
        
        # nodes
        node_header = self.generate_node_header()
        for name, node in self.nodes.iteritems():
            attrs = [None] * len(node_header)
            for i, entry in enumerate(node_header):
                try:
                    attrs[i] = node.info[entry]
                except KeyError:
                    pass
            nodes.append(name+'\t'+'\t'.join(map(str,attrs)))
        with open(self.name+"_nodes.txt", 'w') as f:
            f.write('\n'.join(nodes)+'\n')
        
        # edges
        edge_header = self.generate_edge_header()
        edge_names = self.generate_edges()
        if len(edge_header) > 0:
            edges = list()
            for edge in edge_names:
                attrs = [None] * len(edge_header)
                for i, entry in enumerate(edge_header):
                    try:
                        attrs[i] = self.edge_info[edge][entry]
                    except KeyError:
                        pass
                edges.append('\t'.join(edge)+'\t'+'\t'.join(map(str,attrs)))
        else:
            edges = ['\t'.join(edge_name) for edge_name in edge_names]
            print edges
        
        with open(self.name+"_edges.txt", 'w') as f:
            f.write('\n'.join(edges)+'\n')
    
class Node:
    """class for each node, holds targets and annotation info"""
    def __init__(self, targets=list(), info=dict()):
        self.targets = targets
        self.info = info

    def add_info(self, key, value):
        self.info[key] = value
    
    def add_info_dict(self, info_dict):
        for key, value in info_dict.iteritems():
            add_info(key, value)
    
    def get_degree(self):
        return len(self.targets)