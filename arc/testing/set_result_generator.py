import biopy
import sys
import multiprocessing
import os
import subprocess

DATA_DIR = "/dlab/data/analysis"


def disopred3(fasta_object):
    start = 1
    end = 2
    if len(fasta_object.sequence()) > 18000:
        return "-", "-"
    if not os.path.isdir("%s/%s/DISOPRED3" % (DATA_DIR, fasta_object.organism())):
        os.makedirs("%s/%s/DISOPRED3" % (DATA_DIR, fasta_object.organism()))
    if not os.path.exists(
            "%s/%s/BLAST/%s.mtx" % (DATA_DIR, fasta_object.organism(), fasta_object.accession())):
        # Create the required blast files
        if not os.path.isdir("/dlab/data/analysis/%s" % fasta_object.organism()):
            os.makedirs("/dlab/data/analysis/%s" % fasta_object.organism(), exist_ok=True)
        if not os.path.isdir("/dlab/data/analysis/%s/BLAST" % fasta_object.organism()):
            os.makedirs("/dlab/data/analysis/%s/BLAST" % fasta_object.organism(), exist_ok=True)
        if not os.path.isfile(
                "/dlab/data/analysis/%s/SEQ/%s.fasta" % (fasta_object.organism(), fasta_object.accession())):
            if not os.path.isdir("/dlab/data/analysis/%s/SEQ" % fasta_object.organism()):
                os.makedirs("/dlab/data/analysis/%s/SEQ" % fasta_object.organism(), exist_ok=True)
            f = open("/dlab/data/analysis/%s/SEQ/%s.fasta" % (fasta_object.organism(), fasta_object.accession()), "w+")
            f.write(fasta_object.fasta())
            f.close()

        if not os.path.isfile(
                "/dlab/data/analysis/%s/BLAST/%s.xml" % (
                fasta_object.organism(), fasta_object.accession())) or not os.path.isfile(
            "/dlab/data/analysis/%s/BLAST/%s.chk" % (fasta_object.organism(), fasta_object.accession())):
            proc = subprocess.Popen(
                '/dlab/home/gerdos/bin/blast/bin/blastpgp -a 1 -m 7 -j 3 -h 0.001 -b 5000 -d /dlab/data/UNIPROT/UNIREF90/uniref90.fasta -i /dlab/data/analysis/%s/SEQ/%s.fasta -C /dlab/data/analysis/%s/BLAST/%s.chk' % (
                    fasta_object.organism(), fasta_object.accession(), fasta_object.organism(),
                    fasta_object.accession()), shell=True,
                stdout=subprocess.PIPE)
            result = proc.stdout.read()
            f = open("/dlab/data/analysis/%s/BLAST/%s.xml" % (fasta_object.organism(), fasta_object.accession()), "w")
            f.write(result.decode())
            f.close()
        if not os.path.isfile(
                "/dlab/data/analysis/%s/BLAST/%s.mtx" % (fasta_object.organism(), fasta_object.accession())):
            f = open("/dlab/data/analysis/%s/BLAST/%s.pn" % (fasta_object.organism(), fasta_object.accession()), "w")
            f.write("%s.chk" % fasta_object.accession())
            f.close()
            f = open("/dlab/data/analysis/%s/BLAST/%s.sn" % (fasta_object.organism(), fasta_object.accession()), "w")
            f.write("../SEQ/%s.fasta" % fasta_object.accession())
            f.close()
            proc = subprocess.Popen(
                '/dlab/home/gerdos/bin/blast/bin/makemat -P /dlab/data/analysis/%s/BLAST/%s' % (
                    fasta_object.organism(), fasta_object.accession()), shell=True,
                stdout=subprocess.PIPE)
            proc.communicate()
            os.remove("/dlab/data/analysis/%s/BLAST/%s.pn" % (fasta_object.organism(), fasta_object.accession()))
            os.remove("/dlab/data/analysis/%s/BLAST/%s.sn" % (fasta_object.organism(), fasta_object.accession()))

    disores = []
    transres = []
    try:
        filen = open("%s/%s/DISOPRED3/%s.diso" %
                     (DATA_DIR, fasta_object.organism(), fasta_object.accession()))
        for fil_line in filen:
            if fil_line.startswith("#"):
                continue
            if not fil_line.rstrip():
                continue
            if start <= int(fil_line.split()[0]) <= end:
                disores.append(float(fil_line.split()[3]))
        filen.close()
        filen = open("%s/%s/DISOPRED3/%s.pbdat" %
                     (DATA_DIR, fasta_object.organism(), fasta_object.accession()))
        for fn_line in filen:
            if fn_line.startswith("#"):
                continue
            if not fn_line.rstrip():
                continue
            if start <= int(fn_line.split()[0]) <= end:
                try:
                    transres.append(float(fn_line.split()[3]))
                except ValueError:
                    transres.append(0)

        filen.close()
    except FileNotFoundError:
        proc = subprocess.Popen(
            'perl /dlab/home/gerdos/bin/DISOPRED/run_disopred.pl %s/%s/SEQ/%s.fasta %s/%s/BLAST/%s.mtx 1>/dev/null' % (
                DATA_DIR, fasta_object.organism(), fasta_object.accession(
                ), DATA_DIR, fasta_object.organism(),
                fasta_object.accession()),
            shell=True,
            stdout=subprocess.PIPE)
        if proc.communicate()[1]:  # In case of some error
            return "!", "!"
        try:  # Moving the result files
            os.rename("%s/%s/SEQ/%s.diso" % (DATA_DIR, fasta_object.organism(), fasta_object.accession()),
                      "%s/%s/DISOPRED3/%s.diso" % (
                          DATA_DIR, fasta_object.organism(), fasta_object.accession()))
            os.rename("%s/%s/SEQ/%s.pbdat" % (DATA_DIR, fasta_object.organism(), fasta_object.accession()),
                      "%s/%s/DISOPRED3/%s.pbdat" % (
                          DATA_DIR, fasta_object.organism(), fasta_object.accession()))
        except OSError:
            sys.stderr.write("%s DISOPRED problem\n" %
                             fasta_object.accession())
            return "!", "!"
        filen = open("%s/%s/DISOPRED3/%s.diso" %
                     (DATA_DIR, fasta_object.organism(), fasta_object.accession()))
        for my_line in filen:
            if my_line.startswith("#"):
                continue
            if not my_line.rstrip():
                continue
            if start <= int(my_line.split()[0]) <= end:
                disores.append(float(my_line.split()[3]))
        filen.close()
        filen = open("%s/%s/DISOPRED3/%s.pbdat" %
                     (DATA_DIR, fasta_object.organism(), fasta_object.accession()))
        for my_line in filen:
            if my_line.startswith("#"):
                continue
            if not my_line.rstrip():
                continue
            if start <= int(my_line.split()[0]) <= end:
                try:
                    transres.append(float(my_line.split()[3]))
                except ValueError:
                    transres.append(0)

        filen.close()
    if not transres or not disores:
        try:
            os.unlink("%s/%s/DISOPRED3/%s.diso" % (DATA_DIR, fasta_object.organism(), fasta_object.accession()))
            os.unlink("%s/%s/DISOPRED3/%s.pbdat" % (DATA_DIR, fasta_object.organism(), fasta_object.accession()))
        except FileNotFoundError:
            pass
    print(fasta_object.accession())
    return transres
    # return max(transres)


organism_dct = []
with open(sys.argv[1]) as fn:
    for line in fn:
        if line.startswith('#'):
            continue
        obj = biopy.GetFasta(line.split()[0], server=True)
        # if len(obj.sequence()) > 2000:
        try:
            obj.write()
        except OSError:
            os.unlink('/dlab/data/analysis/{}/SEQ/{}.fasta'.format(obj.organism(), obj.accession()))
            obj.write()
        organism_dct.append(obj)
        # try:
        #     obj.write()
        # except OSError:
        #     os.unlink('/dlab/data/analysis/{}/SEQ/{}.fasta'.format(obj.organism(), obj.fasta()))
        #     obj.write()

print("Reading done!")
p = multiprocessing.Pool(multiprocessing.cpu_count())
for x in p.imap_unordered(disopred3, organism_dct):
    pass
print("Disopred3 done")
p.close()
