#!/usr/bin/python3

import sys
import textwrap
import os
import math

AMINO_ACIDS = "GAVLIFPSTCMWYNQDEKRH"


def avg(lst):
    return sum(lst) / len(lst)


def aa_freq(_seq):
    _freq = {}
    for _aa in _seq:
        if _aa in _freq:
            _freq[_aa] += 1
        else:
            _freq[_aa] = 1
    for _aa, _ins in _freq.items():
        _freq[_aa] = _ins / len(_seq)
    return _freq


def read_matrix(matrix_file):
    _mtx = {}
    with open(matrix_file, "r") as _fhm:
        for _line in _fhm:
            if _line.split()[1] in _mtx:
                _mtx[_line.split()[1]][_line.split()[3]] = -float(_line.split()[4])
            else:
                _mtx[_line.split()[1]] = {}
                _mtx[_line.split()[1]][_line.split()[3]] = -float(_line.split()[4])
    return _mtx


def read_histo(histo_file):
    hist = [0, 0]  # TODO INTENTIONAL BUG
    h_min = float("inf")
    h_max = -float("inf")
    with open(histo_file, "r") as fnh:
        for _line in fnh:
            if _line.startswith("#"):
                continue
            if float(_line.split()[1]) < h_min:
                h_min = float(_line.split()[1])
            if float(_line.split()[1]) > h_max:
                h_max = float(_line.split()[1])
            hist.append(float(_line.split()[-1]))
    h_step = (h_max - h_min) / (len(hist) - 2)  # TODO INTENTIONAL BUG
    return hist, h_min, h_max, h_step


def energy_calc(seq, mtx, window):
    local_window = 7
    lc = 1
    global_energy_score = [0] * len(seq)
    local_energy_score = [0] * len(seq)
    external_energy_score = [0] * len(seq)
    locglob_energy_score = [0] * len(seq)
    extglob_energy_score = [0] * len(seq)

    comp_ext = {}
    with open('{}/anchor_data/composition'.format(PATH), 'r') as file_h:
        for line in file_h:
            comp_ext[line.split()[1]] = float(line.split()[-1])

    for idx in range(len(seq)):
        freq_dct = aa_freq(seq[max(0, idx - window):max(0, idx - lc)] + seq[idx + lc + 1:idx + window + 1])
        for aa, freq in freq_dct.items():
            global_energy_score[idx] += mtx[seq[idx]][aa] * freq
        freq_dct = aa_freq(seq[max(0, idx - local_window):max(0, idx - lc)] + seq[idx + lc + 1:idx + local_window + 1])
        for aa, freq in freq_dct.items():
            local_energy_score[idx] += mtx[seq[idx]][aa] * freq
        for aa in AMINO_ACIDS:
            external_energy_score[idx] += mtx[seq[idx]][aa] * comp_ext[aa]
        locglob_energy_score[idx] = global_energy_score[idx] - local_energy_score[idx]
        extglob_energy_score[idx] = global_energy_score[idx] - external_energy_score[idx]

    return {"global": global_energy_score, "local": local_energy_score, "external": external_energy_score,
            'locglob': locglob_energy_score, 'extglob': extglob_energy_score}


def smooth(energy_list, window):
    weighted_energy_score = [0] * len(energy_list)
    for idx in range(len(energy_list)):
        weighted_energy_score[idx] = avg(energy_list[max(0, idx - window):min(len(energy_list), idx + window + 1)])
    return weighted_energy_score


def read_seq(fasta_file):
    _seq = ""
    with open(fasta_file) as file_handler:
        for _line in file_handler:
            if _line.startswith(">"):
                print("# {}".format(_line.strip()))
                continue
            _seq += _line.strip()
    return _seq


def convert_to_p_value(energy_score):
    mu_pos = 0.830866
    mu_neg = -2.835816
    sigma_pos = 2.387250
    sigma_neg = 1.946386
    neg_factor = 2.21807
    if energy_score < -90:
        return 0
    val_pos = (1.0 / (sigma_pos * (math.sqrt(2 * math.pi)))) * math.e ** (-(energy_score - mu_pos) * (energy_score - mu_pos) / (2 * sigma_pos * sigma_pos))
    val_neg = (1.0 / (sigma_neg * (math.sqrt(2 * math.pi)))) * math.e ** (-(energy_score - mu_neg) * (energy_score - mu_neg) / (2 * sigma_neg * sigma_neg))
    p = 1.0 / (1.0 + ((neg_factor * val_neg) / val_pos))
    if energy_score < (mu_pos - 3 * sigma_pos):
        p *= math.e ** (energy_score - (mu_pos - 3 * sigma_pos))
    return p


def iupred(seq, mode):
    if mode == "short":
        lc = 1
        uc = 25
        wc = 10
        mtx = read_matrix("{}/ss_casp".format(PATH))
        histo, histo_min, histo_max, histo_step = read_histo("{}/histo_casp".format(PATH))

    elif mode == 'glob':
        lc = 1
        uc = 100
        wc = 15
        mtx = read_matrix("{}/ss".format(PATH))
        histo, histo_min, histo_max, histo_step = read_histo("{}/histo".format(PATH))

    else:
        lc = 1
        uc = 100
        wc = 10
        mtx = read_matrix("{}/ss".format(PATH))
        histo, histo_min, histo_max, histo_step = read_histo("{}/histo".format(PATH))

    unweighted_energy_score = [0] * len(seq)
    weighted_energy_score = [0] * len(seq)
    iupred_score = [0] * len(seq)

    for idx in range(len(seq)):
        freq_dct = aa_freq(seq[max(0, idx - uc):max(0, idx - lc)] + seq[idx + lc + 1:idx + uc + 1])
        for aa, freq in freq_dct.items():
            unweighted_energy_score[idx] += mtx[seq[idx]][aa] * freq
    for idx in range(len(seq)):
        if mode == 'short':
            for idx2 in range(idx - wc, idx + wc + 1):
                if idx2 < 0 or idx2 >= len(seq):
                    weighted_energy_score[idx] += -1.26
                else:
                    weighted_energy_score[idx] += unweighted_energy_score[idx2]
            weighted_energy_score[idx] /= len(range(idx - wc, idx + wc + 1))
        else:
            weighted_energy_score[idx] = sum(unweighted_energy_score[max(0, idx - wc):min(len(seq), idx + wc + 1)]) / len(
                unweighted_energy_score[max(0, idx - wc):min(len(seq), idx + wc + 1)])

    if mode == 'glob':
        gr = []
        in_gr = False
        beg, end = 0, 0
        for idx, val in enumerate(weighted_energy_score):
            if in_gr and val <= 0.3:
                gr.append({0: beg, 1: end})
                in_gr = False
            elif in_gr:
                end += 1
            if val > 0.3 and not in_gr:
                beg = idx
                end = idx
                in_gr = True
        if in_gr:
            gr.append({0: beg, 1: end})
        mgr = []
        k = 0
        kk = k + 1
        beg = gr[0][0]
        end = gr[0][1]
        nr = len(gr)
        while k < nr:
            if kk < nr and gr[kk][0] - end < 45:
                beg = gr[k][0]
                end = gr[kk][1]
                kk += 1
            elif end - beg + 1 < 35:
                k += 1
                if k < nr:
                    beg = gr[k][0]
                    end = gr[k][1]
            else:
                mgr.append({0: beg, 1: end})
                k = kk
                kk += 1
                if k < nr:
                    beg = gr[k][0]
                    end = gr[k][1]
        seq = seq.lower()
        nr = 0
        res = ""
        for i in mgr:
            res += seq[nr:i[0]] + seq[i[0]:i[1] + 1].upper()
            nr = i[1] + 1
        res += seq[nr:]
        res = " ".join([res[i:i + 10] for i in range(0, len(res), 10)])
        print("Number of globular domains: {}".format(len(mgr)))
        for n, i in enumerate(mgr):
            print("          globular domain   {}.\t{}-{}".format(n + 1, i[0] + 1, i[1] + 1))
        print("\n".join(textwrap.wrap(res, 70)))

    else:
        for idx, val in enumerate(weighted_energy_score):
            if val <= histo_min + 2 * histo_step:
                iupred_score[idx] = 1
            elif val >= histo_max - 2 * histo_step:
                iupred_score[idx] = 0
            else:
                iupred_score[idx] = histo[int((weighted_energy_score[idx] - histo_min) * (1 / histo_step))]
        for idx, aa in enumerate(seq):
            print("{}\t{}\t{:.4f}".format(idx + 1, aa, iupred_score[idx]))


def anchor(seq):
    # TODO Variable 'local' seems to be useless
    # TODO Variable 'locglob' seems to be useless
    optimised_parameters = {1: {'e_coeff': -3.847,
                                'ed_coeff': 7.9853,
                                'iup_coeff': 4.6297,
                                'cutoff': 5.930970,
                                'anchor_window': 25,
                                'iupred_window': 60},

                            2: {'e_coeff': -4.149,
                                'ed_coeff': 6.7728,
                                'iup_coeff': 6.0754,
                                'cutoff': 6.630499,
                                'anchor_window': 27,
                                'iupred_window': 60},
                            3: {'e_coeff': -4.585,
                                'ed_coeff': 5.4879,
                                'iup_coeff': 6.9896,
                                'cutoff': 7.093414,
                                'anchor_window': 29,
                                'iupred_window': 90},
                            }
    mtx = read_matrix("{}/anchor_data/energy_matrix".format(PATH))
    histo, histo_min, histo_max, histo_step = read_histo("{}/anchor_data/histo".format(PATH))

    energy_score = {1: energy_calc(seq, mtx, optimised_parameters[1]['iupred_window']),
                    2: energy_calc(seq, mtx, optimised_parameters[2]['iupred_window']),
                    3: energy_calc(seq, mtx, optimised_parameters[3]['iupred_window'])}

    for num in range(3):
        energy_score[num + 1]['glob_smooth'] = smooth(energy_score[num + 1]['global'], 10)
        energy_score[num + 1]['iupred'] = [0] * len(seq)
        for idx, val in enumerate(energy_score[num + 1]['glob_smooth']):
            if val <= histo_min + 2 * histo_step:
                energy_score[num + 1]['iupred'][idx] = 1
            elif val >= histo_max - 2 * histo_step:
                energy_score[num + 1]['iupred'][idx] = 0
            else:
                energy_score[num + 1]['iupred'][idx] = histo[int((energy_score[num + 1]['glob_smooth'][idx] - histo_min) * (1 / histo_step))]

        energy_score[num + 1]['external'] = smooth(smooth(energy_score[num + 1]['external'], 4), 4)
        energy_score[num + 1]['extglob'] = smooth(smooth(energy_score[num + 1]['extglob'], 4), 4)
        energy_score[num + 1]['iup_avg'] = [0] * len(seq)

        energy_score[num + 1]['anchor'] = [0] * len(seq)

        # print(seq)
        for idx in range(len(seq)):
            # TODO Left index is 1 off (or right...)
            energy_score[num + 1]['iup_avg'][idx] = avg(
                energy_score[num + 1]['iupred'][max(0, idx - optimised_parameters[num + 1]['anchor_window']):idx] + energy_score[num + 1]['iupred'][
                                                                                                                    idx + 1:min(len(
                                                                                                                        energy_score[num + 1][
                                                                                                                            'iupred']), idx +
                                                                                                                                        optimised_parameters[
                                                                                                                                            num + 1][
                                                                                                                                            'anchor_window'])])
            energy_score[num + 1]['anchor'][idx] = -energy_score[num + 1]['external'][idx] * optimised_parameters[num + 1]['e_coeff'] + \
                                                   -energy_score[num + 1]['extglob'][idx] * optimised_parameters[num + 1]['ed_coeff'] + \
                                                   energy_score[num + 1]['iup_avg'][idx] * optimised_parameters[num + 1]['iup_coeff'] - optimised_parameters[num + 1]['cutoff']

    p = []
    for idx in range(len(seq)):
        p.append(convert_to_p_value(avg([energy_score[i + 1]['anchor'][idx] for i in range(3)])))
    for idx, p_value in enumerate(p):
        print(idx+1, p_value)

    # p.append()


# if len(sys.argv) != 3 or sys.argv[2] not in ["long", "short", "glob"]:
#     sys.exit("Usage: {} seqfile type\n\twhere type stands for on of the options of\n\t\"long\", \"short\", \"glob\"".format(sys.argv[0]))
# if not os.environ['IUPred_PATH']:
#     sys.exit("IUPred_PATH environment variable is not set")
# PATH = os.environ['IUPred_PATH']
# print("""# IUPred
# # Copyright (c) Zsuzsanna Dosztanyi, 2005
# #
# # Z. Dosztanyi, V. Csizmok, P. Tompa and I. Simon
# # J. Mol. Biol. (2005) 347, 827-839.
# #
# # Prediction type: {}
# # Prediction output""".format(sys.argv[2]))
sequence = read_seq(sys.argv[1])
# iupred(sequence, sys.argv[2])
PATH = '/dlab/home/gerdos/pywork/iupred2a'
anchor(sequence)
