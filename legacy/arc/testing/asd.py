import biopy
import random
import re
import os
import sys
import matplotlib.pyplot as plt
import subprocess


def phob(fasta_object, start, end):
    start, end = int(start), int(end)
    # signalp_obj = biopy.SIGNALP(fasta_object)
    # try:
    #     signalp_result = signalp_obj.read()
    # except IOError:
    #     signalp_obj.write()
    #     signalp_result = signalp_obj.read()
    # signal = False
    # for signalp_result_line in signalp_result.splitlines():
    #     if signalp_result_line.startswith("#"):
    #         continue
    #     if not signalp_result_line.rstrip():
    #         continue
    #     if signalp_result_line.split()[9] == "Y":
    #         signal = True

    phobius_obj = biopy.PHOBIUS(fasta_object)
    try:
        phobius_full_result = phobius_obj.read()
    except IOError:
        # return "ETX"
        phobius_obj.write()
        phobius_full_result = phobius_obj.read()
    #
    # if not signal:  # In case of absent signal peptide, its intracellular
    #     return "INT"

    phobius_result = "INT"
    for phobius_result_line in phobius_full_result.splitlines():
        if phobius_result_line.startswith("#"):
            continue
        if not phobius_result_line.rstrip():
            continue
        if re.search("TRANSMEM", phobius_result_line):
            if (int(phobius_result_line.split()[2]) <= start <= int(phobius_result_line.split()[3])) or (
                    int(phobius_result_line.split()[2]) <= end <= int(phobius_result_line.split()[3])):
                phobius_result = "TM"
        if re.search("NON CYTOPLASMIC.", phobius_result_line):
            if (int(phobius_result_line.split()[2]) <= start <= int(phobius_result_line.split()[3])) or (
                    int(phobius_result_line.split()[2]) <= end <= int(phobius_result_line.split()[3])):
                phobius_result = "EXT"

    return phobius_result


def pf(fasta_object):
    pfam_obj = biopy.PFAM(fasta_object)
    try:
        pfam_result = pfam_obj.read()
    except IOError:
        # return "Domain"
        pfam_obj.write()
        pfam_result = pfam_obj.read()
    pfamres = {}
    for pfam_result_line in pfam_result.splitlines():
        if pfam_result_line.startswith("#"):
            continue
        if not pfam_result_line.rstrip():
            continue
        if pfam_result_line.split()[7] == "Domain":
            pfamres[int(pfam_result_line.split()[3])] = int(pfam_result_line.split()[4])

    return pfamres


def get_random_regions_like_dibs():
    human_proteom_acc = []
    with open('/dlab/data/ANCHOR2/human_proteome_reviewed.fasta', "r") as fn:
        for line in fn:
            if line.startswith('>'):
                human_proteom_acc.append(line.split("|")[1])
    with open("/dlab/data/ANCHOR2/data/testing/DIBS-short_testing_set.txt") as fn:
        for line in fn:
            random.shuffle(human_proteom_acc)
            suc = 0
            print("# Similar length random seq for: {}".format(line.strip()))
            for i in human_proteom_acc:
                if suc > 35:
                    break
                try:
                    fasta_obj = biopy.GetFasta(i)
                except ValueError:
                    continue
                try:
                    os.unlink('/dlab/data/analysis/HOMO_SAPIENS/SEQ/{}.fasta'.format(fasta_obj.accession()))
                    fasta_obj.write()
                except FileNotFoundError:
                    fasta_obj.write()
                if len(fasta_obj.sequence()) <= 30:
                    continue
                known_start = (int(line.split()[1].split("-")[0]))
                known_end = (int(line.split()[1].split("-")[1]))
                motif_len = known_end - known_start + 1
                random_start = random.randint(0, len(fasta_obj.sequence()) - motif_len)
                # print(len(fasta_obj.sequence()), random_start)
                for idx in range(random_start, len(fasta_obj.sequence()) - motif_len):
                    if phob(fasta_obj, idx, idx + motif_len) == "INT":
                        if 'Domain' not in pf(fasta_obj, idx, idx + motif_len):
                            print(i, idx + 1, idx + motif_len)
                            if len(fasta_obj.sequence()) < idx + motif_len:
                                sys.exit("WAT")
                            suc += 1
                            break


def anc_new(fasta_object, motif_dct, domain_dct):
    proc = subprocess.Popen(
        'perl /dlab/data/ANCHOR2/dev/ANCHOR_exp/ANCHOR2_exp3.pl /dlab/data/analysis/{}/SEQ/{}.fasta'.format(
            fasta_object.organism(), fasta_object.accession()), shell=True,
        stdout=subprocess.PIPE)
    # Calculate the average of the values
    real = []
    decoy = []
    ancres = proc.communicate()[0]
    for anchor_result_line in ancres.decode().splitlines():
        if anchor_result_line.startswith("#"):
            continue
        if not anchor_result_line.rstrip():
            continue
        trigger = False
        for i, j in domain_dct.items():
            if i <= int(anchor_result_line.split()[0]) <= j:
                trigger = True
        if trigger:
            continue
        for i, j in motif_dct.items():
            if i <= int(anchor_result_line.split()[0]) <= j:
                real.append(float(anchor_result_line.split()[2]))
            else:
                decoy.append(float(anchor_result_line.split()[2]))

    return real, decoy


def anc(fasta_object,  motif_dct, domain_dct):
    anchor_obj = biopy.ANCHOR(fasta_object)
    try:  # If local ANCHOR is availible lets use it
        ancres = anchor_obj.read()
    except IOError:  # In case it isnt, lets create it!
        anchor_obj.write()
        ancres = anchor_obj.read()
    if not ancres:  # In case something is wrong
        sys.stderr.write("{} had some ANCHOR problems".format(
            fasta_object.accession()))
        return "-"
    # Calculate the average of the values
    real = []
    decoy = []
    for anchor_result_line in ancres.splitlines():
        if anchor_result_line.startswith("#"):
            continue
        if not anchor_result_line.rstrip():
            continue
        trigger = False
        for i, j in domain_dct.items():
            if i <= int(anchor_result_line.split()[0]) <= j:
                trigger = True
        if trigger:
            continue
        for i, j in motif_dct.items():
            if i <= int(anchor_result_line.split()[0]) <= j:
                real.append(float(anchor_result_line.split()[2]))
            else:
                decoy.append(float(anchor_result_line.split()[2]))
    return real, decoy


def disopred3(fasta_object, motif_dct, domain_dct):
    if len(fasta_object.sequence()) > 18000:
        return "-", "-"
    if not os.path.isdir("%s/%s/DISOPRED3" % (DATA_DIR, fasta_object.organism())):
        os.makedirs("%s/%s/DISOPRED3" % (DATA_DIR, fasta_object.organism()))
    # if not os.path.exists(
    #         "%s/%s/BLAST/%s.mtx" % (DATA_DIR, fasta_object.organism(), fasta_object.accession())):
    #     # Create the required blast files
    #     biopy.blast(fasta_object.accession())
    real = []
    decoy = []
    try:
        filen = open("%s/%s/DISOPRED3/%s.pbdat" %
                     (DATA_DIR, fasta_object.organism(), fasta_object.accession()))
        for fn_line in filen:
            if fn_line.startswith("#"):
                continue
            if not fn_line.rstrip():
                continue
            trigger = False
            for i, j in domain_dct.items():
                if i <= int(fn_line.split()[0]) <= j:
                    trigger = True
            if trigger:
                continue
            for i, j in motif_dct.items():
                if i <= int(fn_line.split()[0]) <= j:
                    try:
                        real.append(float(fn_line.split()[3]))
                    except ValueError:
                        real.append(0)
                else:
                    try:
                        decoy.append(float(fn_line.split()[3]))
                    except ValueError:
                        decoy.append(0)
        filen.close()
    except FileNotFoundError:
        return [], []
    return real, decoy
    # return max(transres)


def morfchibi(fasta_object, motif_dct, domain_dct):
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
            return morfchibi(fasta_object, motif_dct, domain_dct)
    if not morfres:
        return "-"
    real = []
    decoy = []
    pseudo_line = 1
    for morfchibi_result_line in morfres.splitlines():
        if morfchibi_result_line.startswith("#"):
            continue
        if not morfchibi_result_line.rstrip():
            continue
        trigger = False
        for i, j in domain_dct.items():
            if i <= pseudo_line <= j:
                trigger = True
        if trigger:
            pseudo_line += 1
            continue
        for i, j in motif_dct.items():
            if i <= pseudo_line <= j:
                real.append(round(float(morfchibi_result_line.split()[1]), 2))
            else:
                decoy.append(round(float(morfchibi_result_line.split()[1]), 2))
        pseudo_line += 1
    return real, decoy


def iup(fasta_object, motif_dct, domain_dct):
    iupred_obj = biopy.IUPred(fasta_object)
    try:  # If local IUPred is availible lets use it
        iupred_obj.read()
    except FileNotFoundError:  # In case it isnt, lets create it!
        iupred_obj.write()
    # Calculate the average of the values
    real = []
    decoy = []
    with open('/dlab/data/analysis/{}/IUPRED/{}.dat'.format(fasta_object.organism(), fasta_object.accession()),
              "r") as fn:
        for iupred_result_line in fn:
            if iupred_result_line.startswith("#"):
                continue
            if not iupred_result_line.rstrip():
                continue
            trigger = False
            for i, j in domain_dct.items():
                if i <= int(iupred_result_line.split()[0]) <= j:
                    trigger = True
            if trigger:
                continue
            for i, j in motif_dct.items():
                if i <= int(iupred_result_line.split()[0]) <= j:
                    real.append(float(iupred_result_line.split()[2]))
                else:
                    decoy.append(float(iupred_result_line.split()[2]))
    return real, decoy


def information_gain():
    known_values = {
        'old_anchor': [],
        'disopred3': [],
        'morfchibi': [],
        # 'iupred': [],
        'new_anchor': [],
    }
    decoy_values = {
        'old_anchor': [],
        'disopred3': [],
        'morfchibi': [],
        # 'iupred': [],
        'new_anchor': [],
    }
    n = 1
    # with open("/dlab/data/ANCHOR2/data/testing/ORD-testing_set.txt") as fn:
    motif_dct = {}
    with open("/dlab/data/ANCHOR2/data/testing/DIBS-short_testing_set.txt") as fn:
        # with open("/dlab/home/gerdos/RESULTS/LC8/known_partners.list") as fn:
        for line in fn:
            start = int(line.split()[1].split("-")[0])
            end = int(line.split()[1].split("-")[1])
            if line.split()[0] in motif_dct:
                motif_dct[line.split()[0]][start] = end
            else:
                motif_dct[line.split()[0]] = {}
                motif_dct[line.split()[0]][start] = end
    for i, j in motif_dct.items():
        print('{}\r'.format(n), end="")
        n += 1
        obj = biopy.GetFasta(i)
        domain_dct = pf(obj)
        # domain_dct = {}
        if len(obj.sequence()) > 2000:
            continue
        old_anchor = anc(obj, j, domain_dct)
        new_anch = anc_new(obj, j, domain_dct)
        diso = disopred3(obj, j, domain_dct)
        morf = morfchibi(obj, j, domain_dct)
        # print(i, j)
        # print(domain_dct)
        # print(morf[0])
        # print(morf[1])
        # print(len(morf[0]), len(morf[1]))
        # print(len(obj.sequence()))
        iupred = iup(obj, j, domain_dct)
        [known_values['old_anchor'].append(i) for i in old_anchor[0]]
        [known_values['disopred3'].append(i) for i in diso[0]]
        [known_values['morfchibi'].append(i) for i in morf[0]]
        # [known_values['iupred'].append(i) for i in iupred[0]]
        [known_values['new_anchor'].append(i) for i in new_anch[0]]
        [decoy_values['old_anchor'].append(i) for i in old_anchor[1]]
        [decoy_values['disopred3'].append(i) for i in diso[1]]
        [decoy_values['morfchibi'].append(i) for i in morf[1]]
        # [decoy_values['iupred'].append(i) for i in iupred[1]]
        [decoy_values['new_anchor'].append(i) for i in new_anch[1]]
    return known_values, decoy_values


# def inf_gain(known_values, random_values):
#     concatenated = {}
#     for i in known_values:
#         if i in concatenated:
#             concatenated[i][0] += 1
#         else:
#             concatenated[i] = {0: 1, 1: 0}
#     for i in random_values:
#         if i in concatenated:
#             concatenated[i][1] += 1
#         else:
#             concatenated[i] = {1: 1, 0: 0}
#     conc_sorted = [x[0] for x in sorted(concatenated.items(), key=lambda x: x[0])]
#     conc_book = [x[1] for x in sorted(concatenated.items(), key=lambda x: x[0])]
#     x = []
#     y = []
#     tp, fp = [], []
#     for i in range(len(conc_sorted)):
#         known_not_found = sum([x[0] for x in conc_book[:i]])
#         random_not_found = sum([x[1] for x in conc_book[:i]])
#         # known_not_found /= len(known_values)
#         # random_not_found /= len(random_values)
#         # print(1-known_not_found, random_not_found)
#         infg = biopy.infor_gain([len(random_values) - random_not_found, len(known_values) - known_not_found],
#                                 [random_not_found, known_not_found])
#         if not math.isnan(infg) and infg != float("inf"):
#             y.append(infg)
#             x.append(conc_sorted[i])
#             tp.append((len(known_values) - known_not_found) / len(known_values))
#             fp.append(((len(random_values) - random_not_found) / len(random_values)))
#     print('{:.4f} {:.4f}'.format(tp[y.index(max(y))], fp[y.index(max(y))]))
#     return x, y


DATA_DIR = "/dlab/data/analysis"
# get_random_regions_like_dibs()
# matplotlib.rc('xtick', labelsize=23)
# matplotlib.rc('ytick', labelsize=23)
known_scores, random_scores = information_gain()
color = ['#4477AA','#117733','#CC6677','#88CCEE','#AA4499']
for n, i in enumerate(sorted(known_scores.keys())):
    mrf = biopy.roc(known_scores[i], random_scores[i])
    plt.plot([x[0] for x in sorted(mrf.items())], [x[1] for x in sorted(mrf.items())], label=i, linewidth=2, color=color[n])
plt.plot([0,1], [0,1],'--', color='black')
plt.legend(loc=4)
title = "ROC_dibs_vs_dibs.png"
plt.title(title)
plt.savefig(title, dpi=300)
plt.show()

#
# for i,j in known_scores.items():
#     print('Known', i, len(j))
#     print('Decoy', i, len(random_scores[i]))
#
# print("Reading done")
# print(sum(known_scores['morfchibi']) / len(known_scores['morfchibi']), sum(random_scores['morfchibi'])/len(random_scores['morfchibi']))
#
# anc_old_infg = biopy.information_gain(known_scores['old_anchor'], random_scores['old_anchor'], tp_fp=True)
# plt.plot(anc_old_infg[0], anc_old_infg[1], linewidth=4, color='#4477AA', label="Anchor OLD TP:{:.0f}% FP:{:.0f}%".format(anc_old_infg[2], anc_old_infg[3]))
#
# anc_new_infg = biopy.information_gain(known_scores['new_anchor'], random_scores['new_anchor'], tp_fp=True)
# plt.plot(anc_new_infg[0], anc_new_infg[1], linewidth=4, color='#117733', label="Anchor NEW TP:{:.0f}% FP:{:.0f}%".format(anc_new_infg[2], anc_new_infg[3]))
#
# morfchibi_ingf = biopy.information_gain(known_scores['morfchibi'], random_scores['morfchibi'], tp_fp=True)
# plt.plot(morfchibi_ingf[0], morfchibi_ingf[1], linewidth=4, color='#CC6677', label="Morfchibi TP:{:.0f}% FP:{:.0f}%".format(morfchibi_ingf[2], morfchibi_ingf[3]))
#
# iupred_infg = biopy.information_gain(known_scores['iupred'], random_scores['iupred'], tp_fp=True)
# plt.plot(iupred_infg[0], iupred_infg[1], linewidth=4, color='#88CCEE', label="IUPred TP:{:.0f}% FP:{:.0f}%".format(iupred_infg[2], iupred_infg[3]))
#
# disopred_infg = biopy.information_gain(known_scores['disopred3'], random_scores['disopred3'], tp_fp=True)
# plt.plot(disopred_infg[0], disopred_infg[1], linewidth=4, color='#AA4499', label="DisoPred3 TP:{:.0f}% FP:{:.0f}%".format(disopred_infg[2], disopred_infg[3]))
# title = "DIBS_vs_DIBS_no_pfam_anchor3"
# plt.xlim([-0.02, 1.02])
# plt.ylabel("Information gain", fontsize=14)
# plt.tight_layout()
# plt.title(title)
# # plt.tick_params(
# #     axis='x',  # changes apply to the x-axis
# #     which='both',  # both major and minor ticks are affected
# #     bottom='off',  # ticks along the bottom edge are off
# #     top='off',  # ticks along the top edge are off
# #     labelbottom='off')  # labels along the bottom edge are off
# plt.legend(fontsize=14)
# plt.legend(loc=4)
# plt.savefig("{}.png".format(title), dpi=300)
# plt.show()
