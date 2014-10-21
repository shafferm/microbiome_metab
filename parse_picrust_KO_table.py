"""parse_picrust_KO_table.py"""

FILE_PATH = "ko_13_5_precalculated.tab"

def get_kos():
    return open(FILE_PATH).readline().strip().split()[1:-1]

def get_organism2kos():
    f = open(FILE_PATH)
    kos = f.readline().strip().split()[1:-1]
    org_dict = dict()
    for line in f:
        line = line.strip().split()
        try:
            int(line[0])
        except IndexError:
            print "IndexError"
            print line
            print len(org_dict)
            break
        except ValueError:
            print "ValueError"
            print line[0]
            print len(org_dict)
            break
        otu = line[0]
        otu_kos = set()
        for i, count in enumerate(line[1:-1]):
            if count > 0:
                otu_kos.add(kos[i])
        org_dict[otu] = otu_kos
    return org_dict