"""network.py"""

from collections import defaultdict

class Network:
    """class to hold network information include nodes, edge info and name.
    nodes is a dict of node objects
    edge info keys are tuples used to annotate edges as defined implicitly by the nodes
    TODO: Make nodes into Panada dataframes???
    TODO: change generate_edges and other generates to actual generator functions
    TODO: Make function to read in from cytoscape file to network object
    """
    
    def __init__(self, name, nodes=None, edge_info=None):
        self.name = name
        if nodes == None:
            self.nodes = {}
        else:
            self.nodes = nodes
        if edge_info == None:
            self.edge_info = defaultdict(dict)
        else:
            self.edge_info = edge_info
    
    def merge_network(self, network):
        if network.nodes != None:
            self.merge_node_dict(network.nodes)
        if network.edge_info != None:
            self.add_edge_info_dict(network.edge_info)
            
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

    def add_node(self, name, node):
        self.nodes[name] = node
    
    def add_node_dict(self, node_dict):
        for key, value in node_dict.iteritems():
            self.add_node(key, value)

    def merge_node(self, name, node):
        if name in self.nodes:
            self.nodes[name].add_targets(node.targets)
            for key, value in node.info.iteritems():
                self.nodes[name].add_info(key, value)
        else:
            self.add_node(name, node)
    
    def merge_node_dict(self, node_dict):
        for key, value in node_dict.iteritems():
            self.merge_node(key, value)    

    def add_edge_info(self, edge, name, info):
        self.edge_info[edge][name] = info

    def add_edge_info_dict(self, edge_info_dict):
        for edge, name in edge_info_dict.iteritems():
            for name, info in edge_info_dict[edge].iteritems():
                self.add_edge_info(edge, name, info)
    
    def remove_nodes(self, nodes_to_remove):
        for name, node in self.nodes.iteritems():
            for node_to_remove in nodes_to_remove:
                node.remove_target(node_to_remove)
        for node in nodes_to_remove:
            try:
                del self.nodes[node]
            except:
                pass

    def print_network(self, loc=""):
        edge_names = set()
        nodes = list()

        # nodes
        node_header = list(self.generate_node_header())
        for name, node in self.nodes.iteritems():
            attrs = [None] * len(node_header)
            for i, entry in enumerate(node_header):
                try:
                    attrs[i] = node.info[entry]
                except KeyError:
                    pass
            nodes.append(name+'\t'+'\t'.join(map(str,attrs)))
        with open(loc+self.name+"_nodes.txt", 'w') as f:
            f.write('id\t'+'\t'.join(node_header)+'\n'+'\n'.join(nodes)+'\n')

        # edges
        edge_header = list(self.generate_edge_header())
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

        with open(loc+self.name+"_edges.txt", 'w') as f:
            f.write('source\ttarget\t'+'\t'.join(edge_header)+'\n'+'\n'.join(edges)+'\n')


class Node:
    """class for each node, holds targets and annotation info"""
    def __init__(self, targets=None, info=None):
        if targets == None:
            self.targets = set()
        else:
            self.targets = set(targets)
        if info == None:
            self.info = {}
        else:
            self.info = info

    def add_info(self, key, value):
        self.info[key] = value

    def add_info_dict(self, info_dict):
        for key, value in info_dict.iteritems():
            self.add_info(key, value)

    def add_target(self, target):
        self.targets.add(target)
    
    def add_targets(self, targets):
        self.targets = self.targets | targets

    def remove_target(self, target):
        try:
            self.targets.remove(target)
        except KeyError:
            pass

    def get_out_degree(self):
        return len(self.targets)
