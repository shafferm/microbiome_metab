from collections import defaultdict

def reverse_dict_of_iterables(dict_):
    new = defaultdict(set)
    for key in dict_:
        for value in dict_[key]:
            new[value].add(key)
    return new
    
def make_set_to_set_dict(dict_):
    """Take a dictionary with an interable as a value and a hashable as a key and turn it into a 
    """
    new = defaultdict(set)
    for key in dict_:
        new[frozenset(dict_[key])].add(key)
    return new