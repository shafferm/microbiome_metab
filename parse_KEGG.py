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

def get_reactions():
    """get compounds for each reaction and each KO
    
    """
    f = open("reaction", 'U')
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
        if hasOrtho == True:
            rxn2co[r] = reacts,prods,rev

    return rxn2co

def get_ko2rxns():
    f = open("reaction", 'U')
    f = f.read()
    f = f.strip().split('///')
    ko2rxns = defaultdict(set)
    for entry in f:
        i = 0
        entry = entry.strip().split('\n')
        kos = []
        hasOrtho = False
        
        while i < len(entry):
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
                ko2rxns[ko].add(r)

    return ko2rxns

def get_pathway2kos():
    f = open("reaction", 'U')
    f = f.read()
    f = f.strip().split('///')
    pathway2kos = defaultdict(set)

    for entry in f:
        i = 0
        entry = entry.strip().split('\n')
        kos = []
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
            for path in paths:
                pathway2kos[path] = pathway2kos[path] | set(kos)

    return pathway2kos
    
def get_rxn2kos():
    """"""
    f = open("reaction", 'U')
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
    f = open("reaction", 'U')
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
            
        if hasOrtho == True:
            if len(r) != 6 and r.startswith('R') == False:
                print r
            rxn_names[r] = name

    return rxn_names

def get_ko_names():
    """"""
    f = open("reaction", 'U')
    f = f.read()
    f = f.strip().split('///')
    ko_names = dict()
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

            if new_start == "ORTHOLOGY":
                kos.append(line.strip().split()[0])
                hasOrtho = True
                start = "ORTHOLOGY"
            elif new_start == "":
                if start == "ORTHOLOGY":
                    ko = line.strip().split()[0],' '.join(line.strip().split()[1:])
                    kos.append(ko)
            else:
                start = new_start
            i+=1
        if hasOrtho == True:
            for ko in kos:
                ko_names[ko[0]] = ko[1]
            
    return ko_names
    
def get_co_names():
    """"""
    co_names = dict()
    
    f = open("compound", 'U')
    f = f.read()
    f = f.strip().split('///')
    for entry in f:
        i = 0
        entry = entry.strip().split('\n')

        while i < len(entry):
            new_start = entry[i][:12].strip()
            line = entry[i][12:].strip()

            if new_start == "ENTRY":
                co = line.strip().split()[0]
                name = co
                start = "ENTRY"
            if new_start == "NAME":
                name = line.strip().split(';')[0]
                start = "NAME"
            i+=1
        
        co_names[co] = name
        
    f = open("glycan", 'U')
    f = f.read()
    f = f.strip().split('///')
    for entry in f:
        i = 0
        entry = entry.strip().split('\n')

        while i < len(entry):
            new_start = entry[i][:12].strip()
            line = entry[i][12:].strip()

            if new_start == "ENTRY":
                co = line.strip().split()[0]
                name = co
                start = "ENTRY"
            if new_start == "NAME":
                name = line.strip().split(';')[0]
                start = "NAME"
            i+=1
        
        co_names[co] = name
    
    return co_names

def get_co_counts():
    """get compounds for each reaction and each KO
    
    """
    f = open("reaction", 'U')
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
    
def parse_reaction_mapformula_cos():
    """Creates a CO set from compounds present in reaction_mapformula.lst file.
    """
    f = open("reaction_mapformula.lst", 'U')
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
    f = open("reaction_mapformula.lst", 'U')
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
        
        
    