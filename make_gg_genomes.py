"""make_genome_folder"""

import cPickle
import os

GENOMES_LOC = "gg_genomes/"
OTU_KO_TABLE_LOC = "databases/ko_13_5_precalculated.tab"

def get_genome(otu):
    try:
        genome = cPickle.load(open(GENOMES_LOC+otu, 'rb'))
        return genome
    except:
        print "file does not exist for this otu"
        return None

def main():
    f = open(OTU_KO_TABLE_LOC)
    kos = f.readline().strip().split()[1:-1]
    try:
        os.mkdir(GENOMES_LOC)
    except:
        pass
    for line in f:
        line = line.strip().split()
        try:
            int(line[0])
        except IndexError:
            print "IndexError"
            print line
            continue
        except ValueError:
            print "ValueError"
            print line[0]
            continue
        otu_kos = list()
        for i, count in enumerate(line[1:-1]):
            if float(count) > 0:
                otu_kos.append(kos[i])
        cPickle.dump(otu_kos, open(GENOMES_LOC+line[0], 'wb'))

if __name__ == '__main__':
    main()
