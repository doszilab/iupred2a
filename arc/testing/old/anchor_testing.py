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
        return "ETX"
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
        return "Domain"
        # pfam_obj.write()
        # pfam_result = pfam_obj.read()
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
    with open("/dlab/data/ANCHOR2/data/training/DIBS-short_training_set.txt") as fn:
        for line in fn:
            random.shuffle(human_proteom_acc)
            suc = 0
            print("# Similar length random seq for: {}".format(line.strip()))
            for i in human_proteom_acc:
                if suc > 23:
                    break
                try:
                    fasta_obj = biopy.GetFasta(i)
                except ValueError:
                    continue
                # try:
                #     os.unlink('/dlab/data/analysis/HOMO_SAPIENS/SEQ/{}.fasta'.format(fasta_obj.accession()))
                #     fasta_obj.write()
                # except FileNotFoundError:
                #     fasta_obj.write()
                if len(fasta_obj.sequence()) <= 30:
                    continue
                known_start = (int(line.split()[1].split("-")[0]))
                known_end = (int(line.split()[1].split("-")[1]))
                motif_len = known_end - known_start  +1
                random_start = random.randint(0, len(fasta_obj.sequence())-motif_len)
                # print(len(fasta_obj.sequence()), random_start)
                for idx in range(random_start, len(fasta_obj.sequence())):
                    if phob(fasta_obj, idx, motif_len) == "INT":
                        if 'Domain' not in pf(fasta_obj, idx, motif_len):
                            print(i, idx, idx + motif_len - 1)
                            if len(fasta_obj.sequence()) < motif_len:
                                sys.exit()
                            suc += 1
                            break

#
# def anc_new(fasta_object, protocol, start, end):
#     start, end = int(start), int(end)
#     # print('perl {} 31 20 /dlab/data/analysis/{}/SEQ/{}.fasta'.format(protocol, fasta_object.organism(), fasta_object.accession()))
#     proc = subprocess.Popen(
#         'perl {} 31 20 /dlab/data/analysis/{}/SEQ/{}.fasta'.format(protocol, fasta_object.organism(), fasta_object.accession()), shell=True, stdout=subprocess.PIPE)
#     # Calculate the average of the values
#     ancscore = []
#     ancres = proc.communicate()[0]
#     for anchor_result_line in ancres.decode().splitlines():
#         if anchor_result_line.startswith("#"):
#             continue
#         if not anchor_result_line.rstrip():
#             continue
#         if start <= int(anchor_result_line.split()[0]) <= end:
#             #     if float(anchor_result_line.split()[2]) >= 0.5 and float(anchor_result_line.split()[3]) == 0:
#             #         continue
#             ancscore.append(float(anchor_result_line.split()[2]))
#
#     return sum(ancscore) / (end - start + 1)
#     # return max(ancscore)
#
#
#
# def anc(fasta_object, start, end):
#     start, end = int(start), int(end)
#     anchor_obj = biopy.ANCHOR(fasta_object)
#     try:  # If local ANCHOR is availible lets use it
#         ancres = anchor_obj.read()
#     except IOError:  # In case it isnt, lets create it!
#         anchor_obj.write()
#         ancres = anchor_obj.read()
#     if not ancres:  # In case something is wrong
#         sys.stderr.write("{} had some ANCHOR problems".format(
#             fasta_object.accession()))
#         return "-"
#     # Calculate the average of the values
#     ancscore = []
#     if "Could" in ancres:
#         os.unlink('/dlab/data/analysis/{}/ANCHOR/{}.dat'.format(fasta_object.organism(), fasta_object.accession()))
#         fasta_object.write()
#         return anc(fasta_object, start, end)
#     for anchor_result_line in ancres.splitlines():
#         if anchor_result_line.startswith("#"):
#             continue
#         if not anchor_result_line.rstrip():
#             continue
#         if start <= int(anchor_result_line.split()[0]) <= end:
#             # if float(anchor_result_line.split()[2]) >= 0.5 and float(anchor_result_line.split()[3]) == 0:
#             #     continue
#             ancscore.append(float(anchor_result_line.split()[2]))
#
#     return sum(ancscore) / (end - start + 1)
#     # return max(ancscore)
#
#
# def disopred3(fasta_object, start, end):
#     if len(fasta_object.sequence()) > 18000:
#         return "-", "-"
#     if not os.path.isdir("%s/%s/DISOPRED3" % (DATA_DIR, fasta_object.organism())):
#         os.makedirs("%s/%s/DISOPRED3" % (DATA_DIR, fasta_object.organism()))
#     if not os.path.exists(
#             "%s/%s/BLAST/%s.mtx" % (DATA_DIR, fasta_object.organism(), fasta_object.accession())):
#         # Create the required blast files
#         biopy.blast(fasta_object.accession())
#     disores = []
#     transres = []
#     try:
#         filen = open("%s/%s/DISOPRED3/%s.diso" %
#                      (DATA_DIR, fasta_object.organism(), fasta_object.accession()))
#         for fil_line in filen:
#             if fil_line.startswith("#"):
#                 continue
#             if not fil_line.rstrip():
#                 continue
#             if start <= int(fil_line.split()[0]) <= end:
#                 disores.append(float(fil_line.split()[3]))
#         filen.close()
#         filen = open("%s/%s/DISOPRED3/%s.pbdat" %
#                      (DATA_DIR, fasta_object.organism(), fasta_object.accession()))
#         for fn_line in filen:
#             if fn_line.startswith("#"):
#                 continue
#             if not fn_line.rstrip():
#                 continue
#             if start <= int(fn_line.split()[0]) <= end:
#                 try:
#                     transres.append(float(fn_line.split()[3]))
#                 except ValueError:
#                     transres.append(0)
#
#         filen.close()
#     except FileNotFoundError:
#         proc = subprocess.Popen(
#             'perl /dlab/home/gerdos/bin/DISOPRED/run_disopred.pl %s/%s/SEQ/%s.fasta %s/%s/BLAST/%s.mtx 1>/dev/null' % (
#                 DATA_DIR, fasta_object.organism(), fasta_object.accession(
#                 ), DATA_DIR, fasta_object.organism(),
#                 fasta_object.accession()),
#             shell=True,
#             stdout=subprocess.PIPE)
#         if proc.communicate()[1]:  # In case of some error
#             return "!", "!"
#         try:  # Moving the result files
#             os.rename("%s/%s/SEQ/%s.diso" % (DATA_DIR, fasta_object.organism(), fasta_object.accession()),
#                       "%s/%s/DISOPRED3/%s.diso" % (
#                           DATA_DIR, fasta_object.organism(), fasta_object.accession()))
#             os.rename("%s/%s/SEQ/%s.pbdat" % (DATA_DIR, fasta_object.organism(), fasta_object.accession()),
#                       "%s/%s/DISOPRED3/%s.pbdat" % (
#                           DATA_DIR, fasta_object.organism(), fasta_object.accession()))
#         except OSError:
#             sys.stderr.write("%s DISOPRED problem\n" %
#                              fasta_object.accession())
#             return "!", "!"
#         filen = open("%s/%s/DISOPRED3/%s.diso" %
#                      (DATA_DIR, fasta_object.organism(), fasta_object.accession()))
#         for my_line in filen:
#             if my_line.startswith("#"):
#                 continue
#             if not my_line.rstrip():
#                 continue
#             if start <= int(my_line.split()[0]) <= end:
#                 disores.append(float(my_line.split()[3]))
#         filen.close()
#         filen = open("%s/%s/DISOPRED3/%s.pbdat" %
#                      (DATA_DIR, fasta_object.organism(), fasta_object.accession()))
#         for my_line in filen:
#             if my_line.startswith("#"):
#                 continue
#             if not my_line.rstrip():
#                 continue
#             if start <= int(my_line.split()[0]) <= end:
#                 try:
#                     transres.append(float(my_line.split()[3]))
#                 except ValueError:
#                     transres.append(0)
#
#         filen.close()
#     return sum(transres) / (end - start + 1)
#     # return max(transres)
#
#
# def morfchibi(fasta_object, start, end):
#     start, end = int(start), int(end)
#     morfres = ""
#     # Check the data folder
#     if not os.path.isdir("%s/%s/MORFCHIBI" % (DATA_DIR, fasta_object.organism())):
#         os.makedirs("%s/%s/MORFCHIBI" % (DATA_DIR, fasta_object.organism()))
#     try:
#         flen = open("%s/%s/MORFCHIBI/%s.dat" % (DATA_DIR, fasta_object.organism(), fasta_object.accession()), "r")
#         for morfchibi_result_line in flen:
#             if morfchibi_result_line.startswith("#"):
#                 continue
#             if not morfchibi_result_line.rstrip():
#                 continue
#             morfres += morfchibi_result_line
#         flen.close()
#     except FileNotFoundError:
#         proc = subprocess.Popen(
#             'python3 /dlab/home/gerdos/pywork/SCRIPTS/morfchibi.py %s/%s/SEQ/%s.fasta' % (
#                 DATA_DIR, fasta_object.organism(), fasta_object.accession()),
#             shell=True,
#             stdout=subprocess.PIPE)
#         output = proc.communicate()
#         if output[1]:
#             sys.exit("MOTFCHIBI ERROR: %s" % output[1].decode())
#         else:
#             flen = open("%s/%s/MORFCHIBI/%s.dat" % (DATA_DIR,
#                                                     fasta_object.organism(), fasta_object.accession()), "r")
#             for morfchibi_result_line in flen:
#                 if morfchibi_result_line.startswith("#"):
#                     continue
#                 if not morfchibi_result_line.rstrip():
#                     continue
#                 morfres += morfchibi_result_line
#             flen.close()
#     if not morfres:
#         return "-"
#     morfscore = []
#     pseudo_line = 1  # Morfchibi file structure is horrible as well, so i use a pseudo liner
#     for morfchibi_result_line in morfres.splitlines():
#         if morfchibi_result_line.startswith("#"):
#             continue
#         if not morfchibi_result_line.rstrip():
#             continue
#         if start <= pseudo_line <= end:
#             morfscore.append(float(morfchibi_result_line.split()[1]))
#         pseudo_line += 1
#     return sum(morfscore) / (end - start + 1)
#     # return max(morfscore)
#
#
# def iup(fasta_object, start, end):
#     start, end = int(start), int(end)
#     iupred_obj = biopy.IUPred(fasta_object)
#     try:  # If local IUPred is availible lets use it
#         iupred_obj.read()
#     except FileNotFoundError:  # In case it isnt, lets create it!
#         iupred_obj.write()
#     # Calculate the average of the values
#     iupscore = []
#     with open('/dlab/data/analysis/{}/IUPRED/{}.dat'.format(fasta_object.organism(), fasta_object.accession()), "r") as fn:
#         for iupred_result_line in fn:
#             if iupred_result_line.startswith("#"):
#                 continue
#             if not iupred_result_line.rstrip():
#                 continue
#             if start <= int(iupred_result_line.split()[0]) <= end:
#                 iupscore.append(float(iupred_result_line.split()[2]))
#     return sum(iupscore) / (end - start + 1)
#     # return max(iupscore)
#
#
#
# def information_gain():
#     known_values = {
#         'old_anchor': [],
#         'disopred3': [],
#         'morfchibi': [],
#         'iupred': [],
#         'anchor_linear_single': [],
#         'anchor_sigmoid_single': [],
#     }
#     n = 0
#     with open("/dlab/data/ANCHOR2/data/seq_DIBS-short.txt") as fn:
#         # with open("/dlab/home/gerdos/RESULTS/LC8/known_partners.list") as fn:
#         print("N\tID\t{}\t{}\t{}\t{}\t{}\t{}".format(*(known_values.keys())))
#         for line in fn:
#             n += 1
#             obj = biopy.GetFasta(line.split()[0])
#             start = int(line.split()[-1].split("-")[0])
#             end = int(line.split()[-1].split("-")[1])
#             print('{}\t{}\t'.format(n, obj.accession()), end='')
#             # start = int(line.split()[2])
#             # end = start +8
#             known_values['old_anchor'].append(anc(obj, start, end))
#             known_values['disopred3'].append(disopred3(obj, start, end))
#             known_values['morfchibi'].append(morfchibi(obj, start, end))
#             known_values['iupred'].append(iup(obj, start, end))
#             known_values['anchor_linear_single'].append(anc_new(obj, '/dlab/data/ANCHOR2/dev/ANCHOR_linear_single/devANCHOR_linear_single.pl', start, end))
#             known_values['anchor_sigmoid_single'].append(anc_new(obj, '/dlab/data/ANCHOR2/dev/ANCHOR_sigmoid_single/devANCHOR_sigmoid_single.pl', start, end))
#             print('{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}'.format(*([i[-1] for i in known_values.values()])))
#             # print(n)
#     random_values = {
#         'old_anchor': [],
#         'disopred3': [],
#         'morfchibi': [],
#         'iupred': [],
#         'anchor_linear_single': [],
#         'anchor_sigmoid_single': [],
#     }
#     n = 0
#     print("N\tID\t{}\t{}\t{}\t{}\t{}\t{}".format(*(known_values.keys())))
#     with open("/dlab/home/gerdos/pywork/iupred2a/testing/random_nonglobular_int_segments.list") as fn:
#         # with open("/dlab/home/gerdos/work/results/LC8/inf_gain/random_selected_regions") as fn:
#         for line in fn:
#             n += 1
#             if line.startswith("#") or not line.strip():
#                 continue
#             obj = biopy.GetFasta(line.split()[0])
#             start = int(line.split()[-2])
#             end = int(line.split()[-2]) + int(line.split()[-1])
#             # start = int(line.split()[1])
#             # end = start + 8
#             print('{}\t{}\t'.format(n, obj.accession()), end='')
#             if len(obj.sequence()) > 2000:
#                 continue
#             random_values['old_anchor'].append(anc(obj, start, end))
#             random_values['disopred3'].append(disopred3(obj, start, end))
#             random_values['morfchibi'].append(morfchibi(obj, start, end))
#             random_values['iupred'].append(iup(obj, start, end))
#             random_values['anchor_linear_single'].append(anc_new(obj, '/dlab/data/ANCHOR2/dev/ANCHOR_linear_single/devANCHOR_linear_single.pl', start, end))
#             random_values['anchor_sigmoid_single'].append(anc_new(obj, '/dlab/data/ANCHOR2/dev/ANCHOR_sigmoid_single/devANCHOR_sigmoid_single.pl', start, end))
#             print('{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}'.format(*([i[-1] for i in random_values.values()])))
#     print('Reading done')
#     return known_values, random_values
#
#
# def inf_gain(random_values, known_values):
#     x, y = [], []
#     for i in range(0, int(PREC)):
#         pos = 0
#         neg = 0
#         for val in random_values:
#             if val < i / PREC:
#                 neg += 1
#             else:
#                 pos += 1
#         neg /= len(random_values)
#         pos /= len(random_values)
#         known_pos = 0
#         known_neg = 0
#         for val in known_values:
#             if val < i / PREC:
#                 known_neg += 1
#             else:
#                 known_pos += 1
#         known_neg /= len(known_values)
#         known_pos /= len(known_values)
#         infg = biopy.infor_gain([pos, known_pos], [neg, known_neg])
#         if not math.isnan(infg) and infg != float("inf"):
#             x.append(i / PREC)
#             y.append(infg)
#     return x, y
#
# PREC = 100.0
DATA_DIR = "/dlab/data/analysis"
get_random_regions_like_dibs()
#
# # matplotlib.rc('xtick', labelsize=23)
# # matplotlib.rc('ytick', labelsize=23)
# known_scores, random_scores = information_gain()
#
# plt.plot(*(inf_gain(known_scores['old_anchor'], random_scores['old_anchor'])), linewidth=4, color='#4477AA', label="Anchor OLD (off)")
# plt.plot(*(inf_gain(known_scores['anchor_linear_single'], random_scores['anchor_linear_single'])), linewidth=4, color='#117733', label="Anchor linear single")
# plt.plot(*(inf_gain(known_scores['anchor_sigmoid_single'], random_scores['anchor_sigmoid_single'])), linewidth=4, color='#DDCC77', label="Anchor sigmoid single")
# plt.plot(*(inf_gain(known_scores['morfchibi'], random_scores['morfchibi'])), linewidth=4, color='#CC6677', label="Morfchibi")
# plt.plot(*(inf_gain(known_scores['iupred'], random_scores['iupred'])), linewidth=4, color='#88CCEE', label="IUPred")
# plt.plot(*(inf_gain(known_scores['disopred3'], random_scores['disopred3'])), linewidth=4, color='#AA4499', label="DisoPred3")
#
# title = "DIBs_seq_short_avg"
# plt.xlim([-0.02, 1.02])
# plt.ylabel("Information gain", fontsize=14)
# plt.legend()
# plt.tight_layout()
# plt.title(title)
# # plt.tick_params(
# #     axis='x',  # changes apply to the x-axis
# #     which='both',  # both major and minor ticks are affected
# #     bottom='off',  # ticks along the bottom edge are off
# #     top='off',  # ticks along the top edge are off
# #     labelbottom='off')  # labels along the bottom edge are off
# plt.legend(fontsize=14)
# plt.savefig("{}.png".format(title), dpi=300)
# plt.show()
print(phob(biopy.GetFasta('Q8NG92'), 35, 43))