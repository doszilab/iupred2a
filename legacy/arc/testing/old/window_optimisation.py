import subprocess
import biopy
import multiprocessing
import itertools
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def anc_new(acc, start, end, w1, w2):
    start, end = int(start), int(end)
    proc = subprocess.Popen(
        'perl /dlab/data/ANCHOR2/dev/ANCHOR_single/devANCHOR_single-bounded.pl {} {} 0 0 /dlab/data/analysis/{}/SEQ/{}.fasta {} {}'.format(
            w1, w2, organism_dct[acc], acc, start, end), shell=True,
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
    # return [sum(ancscore) / (end - start + 1)]
    # return max(ancscore)
    return ancscore


def info_gain_opt(w):
    w1, w2 = w
    known_values, random_values = [], []
    with open("/dlab/data/ANCHOR2/data/training/DIBS-short_training_set.txt") as fn:
        for line in fn:
            start = int(line.split()[1].split("-")[0])
            end = int(line.split()[1].split("-")[1])
            # [known_values.append(x) for x in iup(line.split()[0], start, end)]
            [known_values.append(x) for x in anc_new(line.split()[0], start, end, w1, w2)]
    n = 1
    with open("/dlab/home/gerdos/pywork/iupred2a/testing/random_nonglobular_int_segments_old.list") as fn:
        for line in fn:
            if line.startswith("#") or not line.strip():
                continue
            start = int(line.split()[-2])
            end = int(line.split()[-2]) + int(line.split()[-1])
            # [random_values.append(x) for x in iup(line.split()[0], start, end)]
            [random_values.append(x) for x in anc_new(line.split()[0], start, end, w1, w2)]
            if n > 500:
                break
            n += 1
    # Calc the inf gain
    concatenated = {}
    for i in known_values:
        if i in concatenated:
            concatenated[i][0] += 1
        else:
            concatenated[i] = {0: 1, 1: 0}
    for i in random_values:
        if i in concatenated:
            concatenated[i][1] += 1
        else:
            concatenated[i] = {1: 1}
    conc_sorted = [x[0] for x in sorted(concatenated.items(), key=lambda x: x[0])]
    conc_book = [x[1] for x in sorted(concatenated.items(), key=lambda x: x[0])]
    dat = {}
    x = []
    y = []
    for i in range(0, len(conc_sorted)):
        known_not_found = sum([x[0] for x in conc_book[:i] if 0 in x])
        random_not_found = sum([x[1] for x in conc_book[:i] if 1 in x])
        known_not_found /= len(known_values)
        random_not_found /= len(random_values)
        # print(1-known_not_found, random_not_found)
        infg = biopy.infor_gain([1 - random_not_found, 1 - known_not_found], [random_not_found, known_not_found])
        if not math.isnan(infg) and infg != float("inf"):
            dat[infg] = {}
            dat[infg]['known_pos'] = 1 - known_not_found
            dat[infg]['known_neg'] = known_not_found
            dat[infg]['random_pos'] = 1 - random_not_found
            dat[infg]['random_neg'] = random_not_found
            y.append(infg)
            x.append(conc_sorted[i])
            # print('{}\t{}'.format(conc_sorted[i], infg))
    # print(known_values)
    # print(random_values)
    plt.plot(x, y, linewidth=4)
    plt.xlabel("Score")
    plt.ylabel("Information gain")
    plt.title("W1={}, W2={}".format(w1, w2))
    plt.savefig("dlab//home/gerdos/pywork/iupred2a/testing/window_optimisation_0_0/{}_{}.png".format(w1, w2))
    # plt.show()
    print("{} {} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}".format(w1, w2, max(y), x[y.index(max(y))],
                                                                   dat[max(y)]['known_pos'],
                                                                   dat[max(y)]['known_neg'],
                                                                   dat[max(y)]['random_pos'],
                                                                   dat[max(y)]['random_neg']))
    return w1, w2, max(y), x[y.index(max(y))]


organism_dct = {}
with open("/dlab/data/ANCHOR2/data/training/DIBS-short_training_set.txt") as fn:
    for line in fn:
        obj = biopy.GetFasta(line.split()[0], server=True)
        organism_dct[obj.accession()] = obj.organism()
n = 1
with open("/dlab/home/gerdos/pywork/iupred2a/testing/random_nonglobular_int_segments_old.list") as fn:
    for line in fn:
        if line.startswith("#") or not line.strip():
            continue
        obj = biopy.GetFasta(line.split()[0], server=True)
        organism_dct[obj.accession()] = obj.organism()
        if n > 500:
            break
        n += 1
print("Reading done")
window_combinations = list(itertools.product(list(range(5, 51, 2)), list(range(10, 100, 5))))
heatmap = np.zeros(shape=(len(set([x[0] for x in window_combinations])), len(set([x[1] for x in window_combinations]))))

p = multiprocessing.Pool(multiprocessing.cpu_count())
for w1, w2, inf_g, max_place in p.imap_unordered(info_gain_opt, window_combinations):
    heatmap[int((w1 - 5) / 2)][int((w2 - 10) / 5)] = inf_g
p.close()
plt.pcolor(heatmap)
plt.colorbar()
plt.xlim([0, heatmap.shape[1]])
plt.ylim([0, heatmap.shape[0]])
plt.savefig('/dlab/home/gerdos/pywork/iupred2a/testing/window_optimisation_0_0/max_infg_heatmap.png')
