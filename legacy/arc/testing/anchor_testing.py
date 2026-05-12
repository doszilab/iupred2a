import biopy
import random
import re
import os
import sys
import math
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


def pf(fasta_object, start, end):
    start, end = int(start), int(end)
    pfam_obj = biopy.PFAM(fasta_object)
    try:
        pfam_result = pfam_obj.read()
    except IOError:
        # return "Domain"
        pfam_obj.write()
        pfam_result = pfam_obj.read()
    pfamres = "-"
    for pfam_result_line in pfam_result.splitlines():
        if pfam_result_line.startswith("#"):
            continue
        if not pfam_result_line.rstrip():
            continue
        if (int(pfam_result_line.split()[3]) <= start <= int(pfam_result_line.split()[4])) or (
                int(pfam_result_line.split()[2]) <= end <= int(pfam_result_line.split()[4])):
            if pfam_result_line.split()[7] == "Domain" or pfam_result_line.split()[7] == "Family":
                pfamres = pfam_result_line.split()[7]
                pfamres += " "
                pfamres += pfam_result_line.split()[6]

    return "%s" % pfamres


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


def anc_new(fasta_object, start, end):
    start, end = int(start), int(end)
    proc = subprocess.Popen(
        'perl /dlab/data/ANCHOR2/dev/ANCHOR_exp/ANCHOR2_exp3.pl /dlab/data/analysis/{}/SEQ/{}.fasta'.format(
            fasta_object.organism(), fasta_object.accession(), start, end), shell=True,
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


def anc(fasta_object, start, end):
    start, end = int(start), int(end)
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
    ancscore = []
    if "Could" in ancres:
        os.unlink('/dlab/data/analysis/{}/ANCHOR/{}.dat'.format(fasta_object.organism(), fasta_object.accession()))
        fasta_object.write()
        return anc(fasta_object, start, end)
    for anchor_result_line in ancres.splitlines():
        if anchor_result_line.startswith("#"):
            continue
        if not anchor_result_line.rstrip():
            continue
        if start <= int(anchor_result_line.split()[0]) <= end:
            # if float(anchor_result_line.split()[2]) >= 0.5 and float(anchor_result_line.split()[3]) == 0:
            #     continue
            ancscore.append(float(anchor_result_line.split()[2]))
    # return sum(ancscore) / (end - start + 1)
    return ancscore


def disopred3(fasta_object, start, end):
    if len(fasta_object.sequence()) > 18000:
        return "-", "-"
    if not os.path.isdir("%s/%s/DISOPRED3" % (DATA_DIR, fasta_object.organism())):
        os.makedirs("%s/%s/DISOPRED3" % (DATA_DIR, fasta_object.organism()))
    # if not os.path.exists(
    #         "%s/%s/BLAST/%s.mtx" % (DATA_DIR, fasta_object.organism(), fasta_object.accession())):
    #     # Create the required blast files
    #     biopy.blast(fasta_object.accession())
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
    return transres
    # return max(transres)


def morfchibi(fasta_object, start, end):
    start, end = int(start), int(end)
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
            flen = open("%s/%s/MORFCHIBI/%s.dat" % (DATA_DIR,
                                                    fasta_object.organism(), fasta_object.accession()), "r")
            for morfchibi_result_line in flen:
                if morfchibi_result_line.startswith("#"):
                    continue
                if not morfchibi_result_line.rstrip():
                    continue
                morfres += morfchibi_result_line
            flen.close()
    if not morfres:
        return "-"
    morfscore = []
    pseudo_line = 1  # Morfchibi file structure is horrible as well, so i use a pseudo liner
    for morfchibi_result_line in morfres.splitlines():
        if morfchibi_result_line.startswith("#"):
            continue
        if not morfchibi_result_line.rstrip():
            continue
        if start <= pseudo_line <= end:
            morfscore.append(round(float(morfchibi_result_line.split()[1]), 2))
        pseudo_line += 1
    return morfscore


def iup(fasta_object, start, end):
    start, end = int(start), int(end)
    iupred_obj = biopy.IUPred(fasta_object)
    try:  # If local IUPred is availible lets use it
        iupred_obj.read()
    except FileNotFoundError:  # In case it isnt, lets create it!
        iupred_obj.write()
    # Calculate the average of the values
    iupscore = []
    with open('/dlab/data/analysis/{}/IUPRED/{}.dat'.format(fasta_object.organism(), fasta_object.accession()),
              "r") as fn:
        for iupred_result_line in fn:
            if iupred_result_line.startswith("#"):
                continue
            if not iupred_result_line.rstrip():
                continue
            if start <= int(iupred_result_line.split()[0]) <= end:
                iupscore.append(float(iupred_result_line.split()[2]))
    return iupscore
    # return max(iupscore)


def information_gain():
    known_values = {
        'old_anchor': [],
        'disopred3': [],
        'morfchibi': [],
        # 'iupred': [],
        'new_anchor': [],
    }
    n = 1
    # with open("/dlab/data/ANCHOR2/data/testing/ORD-testing_set.txt") as fn:
    with open("/dlab/data/ANCHOR2/data/testing/linkers_for_testing-amalgamated.txt") as fn:
        # with open("/dlab/home/gerdos/RESULTS/LC8/known_partners.list") as fn:
        for line in fn:
            print('{}\r'.format(n), end="")
            n += 1
            obj = biopy.GetFasta(line.split()[0])
            start = int(line.split()[1].split("-")[0])
            end = int(line.split()[1].split("-")[1])
            # start = int(line.split()[2])
            # end = start + 8
            # if len(obj.sequence()) > 2000:
            #     continue
            [known_values['old_anchor'].append(i) for i in anc(obj, start, end)]
            # [known_values['disopred3'].append(i) for i in disopred3(obj, start, end)]
            [known_values['morfchibi'].append(i) for i in morfchibi(obj, start, end)]
            # [known_values['iupred'].append(i) for i in iup(obj, start, end)]
            [known_values['new_anchor'].append(i) for i in anc_new(obj, start, end)]
    random_values = {
        'old_anchor': [],
        'disopred3': [],
        'morfchibi': [],
        # 'iupred': [],
        'new_anchor': [],
    }
    n = 1
    print()
    with open("/dlab/data/ANCHOR2/data/testing/DisProt_complementary_testing_set.txt") as fn:
        # with open("/dlab/home/gerdos/work/results/LC8/inf_gain/random_selected_regions") as fn:
        # with open("/dlab/data/ANCHOR2/data/testing/dibs_like_nonglobular_intracellular_decoy.txt") as fn:
        for line in fn:
            if line.startswith("#") or not line.strip():
                continue
            print('{}\r'.format(n), end="")
            n += 1
            obj = biopy.GetFasta(line.split()[0])
            # if n > 100:
            #     break
            # start = int(line.split()[1])
            # end = int(line.split()[2])
            # start = 0
            # end = int(line.split()[1].split("-")[0])
            # start = int(line.split()[1])
            # end = start + 8
            # if len(obj.sequence()) > 2000:
            #     continue
            start = int(line.split()[1].split("-")[0])
            end = int(line.split()[1].split("-")[1])
            [random_values['old_anchor'].append(i) for i in anc(obj, start, end)]
            # [random_values['disopred3'].append(i) for i in disopred3(obj, start, end)]
            [random_values['morfchibi'].append(i) for i in morfchibi(obj, start, end)]
            # [random_values['iupred'].append(i) for i in iup(obj, start, end)]
            [random_values['new_anchor'].append(i) for i in anc_new(obj, start, end)]
    print()
    return known_values, random_values


DATA_DIR = "/dlab/data/analysis"
# get_random_regions_like_dibs()
# matplotlib.rc('xtick', labelsize=23)
# matplotlib.rc('ytick', labelsize=23)
color = ['#4477AA', '#117733', '#CC6677', '#88CCEE', '#AA4499']
known_scores, random_scores = information_gain()
for n, i in enumerate(sorted(known_scores.keys())):
    mrf = biopy.roc(known_scores[i], random_scores[i])
    plt.plot([x[0] for x in sorted(mrf.items())], [x[1] for x in sorted(mrf.items())], label=i, linewidth=2, color=color[n])
plt.plot([0, 1], [0, 1], '--', color='black')
plt.legend(loc=4)
title = "ROC_dibs_vs_disprot.png"
plt.title(title)
plt.savefig(title, dpi=300)
plt.show()

# for i, j in random_scores.items():
#     random_scores[i] = j*2
# random_scores = random_scores * 2
# for i,j in known_scores.items():
#     print('Known', i, len(j))
#     print('Decoy', i, len(random_scores[i]))
#
# print("Reading done")
# print(sum(known_scores['morfchibi'])/ len(known_scores['morfchibi']), sum(random_scores['morfchibi'])/len(random_scores['morfchibi']))
# anc_old_infg = biopy.information_gain(known_scores['old_anchor'], random_scores['old_anchor'], tp_fp=True)
# plt.plot(anc_old_infg[0], anc_old_infg[1], linewidth=4, color='#4477AA', label="Anchor OLD TP:{:.0f}% FP:{:.0f}%".format(anc_old_infg[2], anc_old_infg[3]))
#
# anc_new_infg = biopy.information_gain(known_scores['new_anchor'], random_scores['new_anchor'], tp_fp=True)
# plt.plot(anc_new_infg[0], anc_new_infg[1], linewidth=4, color='#117733', label="Anchor NEW TP:{:.0f}% FP:{:.0f}%".format(anc_new_infg[2], anc_new_infg[3]))
#
# morfchibi_ingf = biopy.information_gain(known_scores['morfchibi'], random_scores['morfchibi'], tp_fp=True)
# plt.plot(morfchibi_ingf[0], morfchibi_ingf[1], linewidth=4, color='#CC6677', label="Morfchibi TP:{:.0f}% FP:{:.0f}%".format(morfchibi_ingf[2], morfchibi_ingf[3]))
#
# # iupred_infg = biopy.information_gain(known_scores['iupred'], random_scores['iupred'], tp_fp=True)
# # plt.plot(iupred_infg[0], iupred_infg[1], linewidth=4, color='#88CCEE', label="IUPred TP:{:.0f}% FP:{:.0f}%".format(iupred_infg[2], iupred_infg[3]))
#
# disopred_infg = biopy.information_gain(known_scores['disopred3'], random_scores['disopred3'], tp_fp=True)
# plt.plot(disopred_infg[0], disopred_infg[1], linewidth=4, color='#AA4499', label="DisoPred3 TP:{:.0f}% FP:{:.0f}%".format(disopred_infg[2], disopred_infg[3]))
# title = "DIBS_vs_Disprot_anchor_exp3"
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
