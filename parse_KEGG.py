"""parse_KEGG.py
Parses the reaction file from KEGG in order to produce a number of data structures.  This script is
designed so that any of it's functions may be imported and so that any of the data structures it
can generate may be written to pickle's for fast use by other scripts.

pathways are just a 5 digit number, orthology is a 5 digit number preceded by K, reaction is a
5 digit number preceded by an R, compounds are a 5 digit number preceded by a C, and glycans
are a 5 digit number preceded by a G

TODO:  Add main method to make pickles for some/all methods 
"""

from collections import defaultdict
from collections import namedtuple

COMPOUND_LOC = "databases/compound"
GLYCAN_LOC = "databases/glycan"
REACTION_LOC = "databases/reaction"
KO_LOC = "databases/ko"
REACTION_MAPFORMULA_LOC = "databases/reaction_mapformula.lst"

class KEGG_Parser:
    """class to retrieve from dictionaries in KEGG"""
    def __init__(self):
        self.rxn2cos = None
        self.ko2rxns = None
        self.pathway2kos = None
        self.pathway2rxns = None
        self.rxn2kos = None
        self.rxn_names = None
        self.ko_names = None
        self.co_names = None
        self.rxns = None
    
    def get_rxns_from_ko(self, ko):
        if self.ko2rxns == None:
            self.ko2rxns = get_ko2rxns()
        try:
            return self.ko2rxns[ko]
        except KeyError:
            print "KO id " + ko + " doesn't exist in this set."
            return set()
    
    def get_co_info(self, co):
        if self.co_names == None:
            self.co_names = get_co_info()
        try:
            return self.co_names[co]
        except KeyError:
            print "CO id " + co + " doesn't exist in this set"
            return None
    
    def get_kos_from_pathway(self, pathway):
        if self.pathway2kos == None:
            self.pathway2kos = get_pathway2kos()
        try:
            return self.pathway2kos[pathway[-5:]]
        except KeyError:
            print "pathway number " + pathway[-5:] + " doesn't exist in this set"
            return set()
            
    def get_rxns_from_pathway(self, pathway):
        if self.pathway2rxns == None:
            self.pathway2rxns = get_pathway2rxns()
        try:
            return self.pathway2rxns[pathway[-5:]]
        except KeyError:
            print "pathway number " + pathway[-5:] + " doesn't exist in this set"
            return set()
            
    def get_kos_from_rxn(self, rxn):
        if self.rxn2kos == None:
            self.rxn2kos = get_rxn2kos()
        try:
            return self.rxn2kos[rxn]
        except KeyError:
            print "reaction id " + rxn + " doesn't exist in this set"
            return set()
    
    def get_rxn(self, rxn):
        if self.rxn2cos == None:
            self.rxn2cos = get_reactions()
        try:
            return self.rxn2cos[rxn]
        except KeyError:
            print "reaction id " + rxn + " doesn't exist in this set"
            return None
    
    def get_rxn_name(self, rxn):
        if self.rxn_names == None:
            self.rxn_names = get_rxn_names()
        try:
            return self.rxn_names[rxn]
        except KeyError:
            print "reaction id " + rxn + " doesn't exist in this set"
            return None

    def get_ko_name(self, ko):
        if self.ko_names == None:
            self.ko_names = get_ko_names()
        try:
            return self.ko_names[ko]
        except KeyError:
            print "KO id " + ko + " doesn't exist in this set."
            return set()

def get_reactions():
    """get compounds for each reaction and each KO
    
    """
    f = open(REACTION_LOC, 'U')
    f = f.read()
    f = f.strip().split('///')
    rxn2co = dict()
    for entry in f:
        i = 0
        entry = entry.strip().split('\n')
        kos = []
        paths = []
        hasOrtho = False
        rev = False
        
        while i < len(entry):
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
            elif new_start == "ORTHOLOGY":
                kos.append(line.strip().split()[0])
                hasOrtho = True
                start = "ORTHOLOGY"
            elif new_start == "":
                if start == "ORTHOLOGY":
                    kos.append(line.strip().split()[0])
            else:
                start = new_start
            i+=1
        rxn2co[r] = reacts,prods,rev

    return rxn2co

# def get_ko2rxns():
#     f = open(REACTION_LOC, 'U')
#     f = f.read()
#     f = f.strip().split('///')
#     ko2rxns = defaultdict(set)
#     for entry in f:
#         i = 0
#         entry = entry.strip().split('\n')
#         kos = []
#         hasOrtho = False
#
#         while i < len(entry):
#             new_start = entry[i][:12].strip()
#             line = entry[i][12:].strip()
#
#             if new_start == "ENTRY":
#                 r = line.strip().split()[0]
#                 start = "ENTRY"
#             elif new_start == "ORTHOLOGY":
#                 kos.append(line.strip().split()[0])
#                 hasOrtho = True
#                 start = "ORTHOLOGY"
#             elif new_start == "":
#                 if start == "ORTHOLOGY":
#                     kos.append(line.strip().split()[0])
#             else:
#                 start = new_start
#             i+=1
#         if hasOrtho == True:
#             if len(r) != 6 and r.startswith('R') == False:
#                 print r
#             for ko in kos:
#                 ko2rxns[ko].add(r)
#
#     return ko2rxns

def get_ko2rxns():
    f = open(KO_LOC, 'U')
    f = f.read()
    f = f.strip().split('///')
    ko2rxns = {}
    for entry in f:
        i = 0
        entry = entry.strip().split('\n')
        ko = None
        rxns = None
        hasOrtho = False
        
        while i < len(entry):
            new_start = entry[i][:12].strip()
            line = entry[i][12:].strip()

            if new_start == "ENTRY":
                ko = line.strip().split()[0]
                start = "ENTRY"
            elif new_start == "DBLINKS":
                line = line.split()
                if line[0] == "RN:":
                    rxns = line[1:]
            elif new_start == "":
                if start == "DBLINKS":
                    if line[0] == "RN:":
                        rxns = line[1:]
            else:
                start = new_start
            i+=1
        if rxns != None:
            ko2rxns[ko] = set(rxns)
        else:
            ko2rxns[ko] = set()

    return ko2rxns

# def get_pathway2kos():
#     f = open(REACTION_LOC, 'U')
#     f = f.read()
#     f = f.strip().split('///')
#     pathway2kos = defaultdict(set)
#
#     for entry in f:
#         i = 0
#         entry = entry.strip().split('\n')
#         kos = []
#         paths = []
#         hasOrtho = False
#
#         while i < len(entry):
#             new_start = entry[i][:12].strip()
#             line = entry[i][12:].strip()
#
#             if new_start == "ENTRY":
#                 r = line.strip().split()[0]
#                 start = "ENTRY"
#             elif new_start == "PATHWAY":
#                 paths.append(line.strip().split()[0][-5:])
#                 start = "PATHWAY"
#             elif new_start == "ORTHOLOGY":
#                 kos.append(line.strip().split()[0])
#                 hasOrtho = True
#                 start = "ORTHOLOGY"
#             elif new_start == "":
#                 if start == "ORTHOLOGY":
#                     kos.append(line.strip().split()[0])
#                 if start == "PATHWAY":
#                     paths.append(line.strip().split()[0][-5:])
#             else:
#                 start = new_start
#             i+=1
#         if hasOrtho == True:
#             if len(r) != 6 and r.startswith('R') == False:
#                 print r
#             for path in paths:
#                 pathway2kos[path] = pathway2kos[path] | set(kos)
#
#     return pathway2kos

def get_pathway2kos():
    f = open(KO_LOC, 'U')
    f = f.read()
    f = f.strip().split('///')
    pathway2kos = defaultdict(set)

    for entry in f:
        i = 0
        entry = entry.strip().split('\n')
        paths = []
        hasOrtho = False
        
        while i < len(entry):
            new_start = entry[i][:12].strip()
            line = entry[i][12:].strip()

            if new_start == "ENTRY":
                ko = line.strip().split()[0]
                start = "ENTRY"
            elif new_start == "PATHWAY":
                paths.append(line.strip().split()[0][-5:])
                start = "PATHWAY"
            elif new_start == "":
                if start == "PATHWAY":
                    paths.append(line.strip().split()[0][-5:])
            else:
                start = new_start
            i+=1
        for path in paths:
            pathway2kos[path].add(ko)

    return pathway2kos

def get_pathway2rxns():
    f = open(REACTION_LOC, 'U')
    f = f.read()
    f = f.strip().split('///')
    pathway2rxns = defaultdict(set)

    for entry in f:
        i = 0
        entry = entry.strip().split('\n')
        paths = []
        hasOrtho = False

        while i < len(entry):
            new_start = entry[i][:12].strip()
            line = entry[i][12:].strip()

            if new_start == "ENTRY":
                r = line.strip().split()[0]
                start = "ENTRY"
            elif new_start == "PATHWAY":
                paths.append(line.strip().split()[0][-5:])
                start = "PATHWAY"
            elif new_start == "":
                if start == "PATHWAY":
                    paths.append(line.strip().split()[0][-5:])
            else:
                start = new_start
            i+=1
        for path in paths:
            pathway2rxns[path].add(r)

    return pathway2rxns
    
def get_rxn2kos():
    """"""
    f = open(REACTION_LOC, 'U')
    f = f.read()
    f = f.strip().split('///')
    rxn2kos = defaultdict(set)
    for entry in f:
        i = 0
        entry = entry.strip().split('\n')
        kos = []
        hasOrtho = False
        
        while i < len(entry):
            rev = False
            new_start = entry[i][:12].strip()
            line = entry[i][12:].strip()

            if new_start == "ENTRY":
                r = line.strip().split()[0]
                start = "ENTRY"
            elif new_start == "ORTHOLOGY":
                kos.append(line.strip().split()[0])
                hasOrtho = True
                start = "ORTHOLOGY"
            elif new_start == "":
                if start == "ORTHOLOGY":
                    kos.append(line.strip().split()[0])
            else:
                start = new_start
            i+=1
            
        if hasOrtho == True:
            if len(r) != 6 and r.startswith('R') == False:
                print r
            for ko in kos:
                rxn2kos[r].add(ko)

    return rxn2kos
    
def get_rxn_names():
    """"""
    f = open(REACTION_LOC, 'U')
    f = f.read()
    f = f.strip().split('///')
    rxn_names = dict()
    for entry in f:
        i = 0
        entry = entry.strip().split('\n')
        kos = []
        hasOrtho = False
        
        while i < len(entry):
            rev = False
            new_start = entry[i][:12].strip()
            line = entry[i][12:].strip()

            if new_start == "ENTRY":
                r = line.strip().split()[0]
                name = r
                start = "ENTRY"
            elif new_start == "NAME":
                name = line.strip()
            elif new_start == "ORTHOLOGY":
                hasOrtho = True
                start = "ORTHOLOGY"
            else:
                start = new_start
            i+=1
            
        if len(r) != 6 and r.startswith('R') == False:
            print r
        rxn_names[r] = name

    return rxn_names

def get_ko_names():
    """"""
    f = open(KO_LOC, 'U')
    f = f.read()
    f = f.strip().split('///')
    ko_names = dict()
    for entry in f:
        i = 0
        entry = entry.strip().split('\n')
        ko = None
        name = None
        defin = None

        while i < len(entry):
            new_start = entry[i][:12].strip()
            line = entry[i][12:].strip()
            
            if new_start == "ENTRY":
                ko = line.strip().split()[0]
            if new_start == "NAME":
                name = line.strip()
            if new_start == "DEFINITION":
                defin = " ".join(line.split()[:-1])
            i+=1
        
        if name == None:
            name = defin
        elif defin != None and name != None:
            name = defin + " (" + name + ")"
        elif defin == None and name == None:
            name = ko
                
        ko_names[ko] = name
            
    return ko_names
    
def get_co_info():
    """returns a named tuple with all avaliable compound data
    nametuple fields: id, name, formula, mass
    """
    co_names = dict()
    Compound = namedtuple('Compound', 'co, name, formula, mass, rxn')
    
    f = open(COMPOUND_LOC, 'U')
    f = f.read()
    f = f.strip().split('///')
    for entry in f:
        i = 0
        entry = entry.strip().split('\n')
        co = None
        name = None
        formula = None
        mass = None

        while i < len(entry):
            new_start = entry[i][:12].strip()
            line = entry[i][12:].strip()

            if new_start == "ENTRY":
                co = line.strip().split()[0]
            if new_start == "NAME":
                name = line.strip().split(';')[0]
            if new_start == "FORMULA":
                formula = line.strip()
            if new_start == "EXACT_MASS":
                mass = line.strip()
            if new_start == "REACTION":
                rxns = line.strip().split()
            i+=1
        
        co_names[co] = Compound(co = co, name = name, formula = formula, mass = mass, rxn = rxns)
        
    f = open(GLYCAN_LOC, 'U')
    f = f.read()
    f = f.strip().split('///')
    for entry in f:
        i = 0
        entry = entry.strip().split('\n')
        co = None
        name = None
        formula = None
        mass = None

        while i < len(entry):
            new_start = entry[i][:12].strip()
            line = entry[i][12:].strip()

            if new_start == "ENTRY":
                co = line.strip().split()[0]
            if new_start == "NAME":
                name = line.strip().split(';')[0]
            if new_start == "COMPOSITION":
                formula = line.strip()
            if new_start == "MASS":
                mass = line.strip().split()[0]
            if new_start == "REACTION":
                rxns = line.strip().split()
            i+=1
        
        co_names[co] = Compound(co = co, name = name, formula = formula, mass = mass, rxn = rxns)
    
    return co_names

def get_co_counts():
    """get compounds for each reaction and each KO
    
    """
    f = open(REACTION_LOC, 'U')
    f = f.read()
    f = f.strip().split('///')
    co_counts = Counter()
    for entry in f:
        i = 0
        entry = entry.strip().split('\n')
        
        while i < len(entry):
            new_start = entry[i][:12].strip()
            line = entry[i][12:].strip()
            
            if new_start == "EQUATION":
                equ = line.split('=>')
                if equ[0][-1] == '<':
                    equ[0] = equ[0][:-1]
                for part in equ[0].strip().split():
                    if part[0] == 'C' or part[0] == 'G':
                        co_counts[part[:6]]+=1
                for part in equ[1].strip().split():
                    if part[0] == 'C' or part[0] == 'G':
                        co_counts[part[:6]]+=1
            i+=1
    return co_counts
    
def get_reaction_mapformula_cos():
    """Creates a CO set from compounds present in reaction_mapformula.lst file.
    """
    f = open(REACTION_MAPFORMULA_LOC, 'U')
    cos = set()
    for line in f:
        line = line.strip().split()[2:]
        for part in line:
            if len(part) == 6 and part.startswith("C"):
                cos.add(part)
    return cos

def parse_reaction_mapformula():
    """adapted from parse_formula() from run_metabolic_networks_old.py
    """
    f = open(REACTION_MAPFORMULA_LOC, 'U')
    rxns = dict()
    for line in f:
        #from parse_mapformula_file from parse_kegg.py
        line = line.strip().split(':')
        rxn = line[0].strip()
        formula = line[2].strip()
        
        
        formula = formula.split('=')
        #get compounds on left and right side of the equation
        left = formula[0][:-1]
        left = left.split('+')
        left = [i.strip() for i in left]
        right = formula[1][1:]
        right = right.split('+')
        right = [i.strip() for i in right]
        rev = False
        
        #Determine whether reaction goes in forward and/or reverse direction
        # and create node objects
        if formula[0].endswith('<'):
            rev = True
        
        rxns[rxn] = left, right, rev
        
    return rxns
