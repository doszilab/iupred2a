import biopy
import subprocess
import sys
import matplotlib.pyplot as plt
import os

DATA_DIR = "/dlab/data/analysis"


def disopred3(fasta_object):
    if len(fasta_object.sequence()) > 18000:
        return "-", "-"
    if not os.path.isdir("%s/%s/DISOPRED3" % (DATA_DIR, fasta_object.organism())):
        os.makedirs("%s/%s/DISOPRED3" % (DATA_DIR, fasta_object.organism()))
    # if not os.path.exists(
    #         "%s/%s/BLAST/%s.mtx" % (DATA_DIR, fasta_object.organism(), fasta_object.accession())):
    #     # Create the required blast files
    #     biopy.blast(fasta_object.accession())
    real = []
    try:
        filen = open("%s/%s/DISOPRED3/%s.pbdat" %
                     (DATA_DIR, fasta_object.organism(), fasta_object.accession()))
        for fn_line in filen:
            if fn_line.startswith("#"):
                continue
            if not fn_line.rstrip():
                continue
            try:
                real.append(float(fn_line.split()[3]))
            except ValueError:
                real.append(0)
        filen.close()
    except FileNotFoundError:
        return []
    return real
    # return max(transres)


def morfchibi(fasta_object):
    morfres = ""
    # Check the data folder
    if not os.path.isdir("%s/%s/MORFCHIBI" % (DATA_DIR, fasta_object.organism())):
        os.makedirs("%s/%s/MORFCHIBI" % (DATA_DIR, fasta_object.organism()))
    try:
        flen = open("%s/%s/MORFCHIBI/%s.dat" % (DATA_DIR, fasta_object.organism(), fasta_object.accession()), "r")
        for morfchibi_result_line in flen:
            if morfchibi_result_line.startswith("#"):
                continue
            if not morfchibi_result_line.rstrip():
                continue
            morfres += morfchibi_result_line
        flen.close()
    except FileNotFoundError:
        proc = subprocess.Popen(
            'python3 /dlab/home/gerdos/pywork/SCRIPTS/morfchibi.py %s/%s/SEQ/%s.fasta' % (
                DATA_DIR, fasta_object.organism(), fasta_object.accession()),
            shell=True,
            stdout=subprocess.PIPE)
        output = proc.communicate()
        if output[1]:
            sys.exit("MOTFCHIBI ERROR: %s" % output[1].decode())
        else:
            return morfchibi(fasta_object)
    if not morfres:
        return "-"
    real = []
    pseudo_line = 1
    for morfchibi_result_line in morfres.splitlines():
        if morfchibi_result_line.startswith("#"):
            continue
        if not morfchibi_result_line.rstrip():
            continue
        real.append(round(float(morfchibi_result_line.split()[1]), 2))
        pseudo_line += 1
    return real


def anc_new(fasta_object):
    # start, end = int(start), int(end)
    proc = subprocess.Popen(
        'perl /dlab/data/ANCHOR2/dev/ANCHOR_single/ANCHOR2_raw.pl /dlab/data/analysis/{}/SEQ/{}.fasta'.format(
            fasta_object.organism(), fasta_object.accession()), shell=True,
        stdout=subprocess.PIPE)
    # Calculate the average of the values
    ancscore = []
    ancres = proc.communicate()[0]
    for anchor_result_line in ancres.decode().splitlines():
        if anchor_result_line.startswith("#"):
            continue
        if not anchor_result_line.rstrip():
            continue
        # if start <= int(anchor_result_line.split()[0]) <= end:
        ancscore.append(float(anchor_result_line.split()[2]))
    return ancscore


obj = biopy.GetFasta(sys.argv[1])
ancy = anc_new(obj)
morfy = morfchibi(obj)
disoy = disopred3(obj)
plt.title(sys.argv[1])
plt.plot(ancy, linewidth=4, label='Anchor')
plt.plot(morfy, linewidth=4, label='Morfchibi')
plt.plot(disoy, linewidth=4, label='Disopred3')
plt.plot([0,len(obj.sequence())], [0.5,0.5], color='black')
plt.legend(loc=4)
# plt.xlim([1100, 1150])
plt.show()


