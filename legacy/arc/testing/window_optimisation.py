import subprocess
import biopy
import multiprocessing




def anc_new(fasta_object, start, end, w):
    start, end = int(start), int(end)
    proc = subprocess.Popen(
        'perl /home/gerdos/ANCHOR2_final-norm.pl /dlab/data/analysis/{}/SEQ/{}.fasta {}'.format(
            fasta_object.organism(), fasta_object.accession(), w), shell=True,
        stdout=subprocess.PIPE)
    # Calculate the average of the values
    ancscore = []
    ancres = proc.communicate()[0]
    for anchor_result_line in ancres.decode().splitlines():
        if anchor_result_line.startswith("#"):
            continue
        if not anchor_result_line.rstrip():
            continue
        if start <= int(anchor_result_line.split()[0]) <= end:
            ancscore.append(float(anchor_result_line.split()[2]))
    return ancscore


def mp(i):
    real_vals = []
    decoy_values = []
    with open("/dlab/data/ANCHOR2/data/training/DIBS-short_training_set.txt") as fn:
        for line in fn:
            obj = biopy.GetFasta(line.split()[0], server=False)
            start = int(line.split()[1].split("-")[0])
            end = int(line.split()[1].split("-")[1])
            real_vals += anc_new(obj, start, end, i)
            # obj.write()
    with open("/dlab/data/ANCHOR2/data/training/decoy_non-glob_int_training.txt") as fn:
        for line in fn:
            if line.startswith("#") or not line.strip():
                continue
            obj = biopy.GetFasta(line.split()[0], server=False)
            start = int(line.split()[1])
            end = int(line.split()[2])
            decoy_values += anc_new(obj, start, end, i)
    print(i, max(biopy.information_gain(real_vals, decoy_values)[1]))

p = multiprocessing.Pool(multiprocessing.cpu_count())
list(p.imap_unordered(mp, list(range(40, 51))))



