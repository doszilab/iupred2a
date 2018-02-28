import biopy
import sys
import os

with open(sys.argv[1]) as fn:
    for line in fn:
        if line.startswith('#'):
            continue
        obj = biopy.GetFasta(line.split()[0], server=False)
        os.system('ln -s /dlab/data/analysis/{0}/DISOPRED3/{1}.pbdat /dlab/data/ANCHOR2/DISOPRED3/{1}.pbdat'.format(obj.organism(), obj.accession()))
        os.system('ln -s /dlab/data/analysis/{0}/DISOPRED3/{1}.diso /dlab/data/ANCHOR2/DISOPRED3/{1}.dat'.format(obj.organism(), obj.accession()))
        os.system('ln -s /dlab/data/analysis/{0}/MORFCHIBI/{1}.dat /dlab/data/ANCHOR2/MORFCHIBI/{1}.dat'.format(obj.organism(), obj.accession()))