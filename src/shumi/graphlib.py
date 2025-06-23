#  Copyright (c) 2025 National Institutes of Health
#  Written by Pierce Radecki
#  This program comes with ABSOLUTELY NO WARRANTY; it is intended for
#  Research Use Only and not for use in diagnostic procedures.

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator

from .pcrlib import MutationEvent, RecombinationEvent

matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 6})

logo = r"""
//////////////////////////////////////////
//                                      //
//           __     __  __ __  ___ ____ //
//    _____ / /_   / / / //  |/  //  _/ //
//   / ___// __ \ / / / // /|_/ / / /   //
//  (__  )/ / / // /_/ // /  / /_/ /    //
// /____//_/ /_/ \____//_/  /_//___/    //
//                                      //
//////////////////////////////////////////
"""


def summary_graph(counts, samples, samples_ccs, fp_out):
    counts = np.array(counts)
    samples = np.array(samples)
    props = counts / counts.sum()
    props_sampled = samples / samples.sum()

    plt.figure(figsize=(12, 10))

    plt.subplot(331)
    if len(set(counts)) == 1:
        # Stem plot if single cDNA
        plt.plot([counts[0], counts[0]], [0, len(counts)], ls='-', color='tomato')
        plt.scatter([counts[0]], [len(counts)], marker='o', s=30, alpha=1, color='tomato')
    else:
        plt.hist(counts, rwidth=2, color='tomato')
    plt.xlabel('Number of PCR DNA molecules generated from bin')
    plt.gca().set_xlim(left=0)
    plt.gca().set_ylim(bottom=0)
    plt.ylabel('Counts')

    plt.subplot(332)
    plt.scatter(counts, samples, marker='o', s=30, alpha=0.5, color='cadetblue')
    plt.plot([0, max(counts)], [0, max(props) * sum(samples)], ls='--', color='k')
    plt.xlabel('Number of PCR DNA molecules generated from bin')
    plt.gca().set_xlim(left=0)
    plt.ylabel('Number of PCR DNA molecules sampled from bin')

    plt.subplot(335)
    plt.hist(samples, bins=50, color='mediumorchid')
    plt.xlabel('Number of PCR DNA molecules sampled from bin')
    plt.gca().set_xlim(left=0)
    plt.ylabel('Counts')

    if samples_ccs:
        samples_ccs = np.array(samples_ccs)
        plt.subplot(333)
        plt.scatter(counts, samples_ccs, marker='o', s=30, alpha=0.5, color='cadetblue')
        plt.plot([0, max(counts)], [0, max(props) * sum(samples_ccs)], ls='--', color='k')
        plt.xlabel('Number of PCR DNA molecules generated from bin')
        plt.gca().set_xlim(left=0)
        plt.ylabel('Number of CCS obtained from bin')

        plt.subplot(336)
        plt.scatter(samples, samples_ccs, marker='o', s=30, alpha=0.5, color='cadetblue')
        plt.plot([0, max(samples)], [0, max(props_sampled) * sum(samples_ccs)], ls='--', color='k')
        plt.xlabel('Number of PCR DNA molecules sampled from bin')
        plt.gca().set_xlim(left=0)
        plt.ylabel('Number of CCS obtained from bin')

        plt.subplot(339)
        plt.hist(samples_ccs, bins=50, color='lightseagreen')
        plt.xlabel('Number of CCS obtained from bin')
        # plt.gca().set_xlim(left=0)
        plt.ylabel('Counts')

    else:
        pass

    plt.tight_layout()
    plt.savefig(fp_out, format='png')


def pcr_summary_graph(simulator, fp_out):
    num_cycles = simulator.num_cycles
    plt.figure(figsize=(10, 10))

    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))

    cops = [mol.times_copied for mol in simulator.current_pool.values()]
    copies = np.bincount(cops)

    muts = [sum(1 if type(e) is MutationEvent else 0 for e in mol.mutation_history) for mol in
            simulator.current_pool.values()]

    mut_counts = np.bincount([sum(1 if type(e) is MutationEvent else 0 for e in mol.mutation_history) for mol in
                              simulator.current_pool.values()])
    recombs = [sum(1 if type(e) is RecombinationEvent else 0 for e in mol.mutation_history) for mol in
               simulator.current_pool.values()]

    rec_counts = np.bincount(recombs)
    plt.subplot(221)
    plt.bar(range(len(copies)), copies)
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.xlabel('Number of times copied')
    plt.ylabel('Molecules')
    plt.xlim((-1, num_cycles + 1))

    plt.subplot(222)
    plt.bar(range(len(mut_counts)), mut_counts)
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.xlabel('Number of mutations')
    plt.ylabel('Molecules')
    plt.xlim((-1, num_cycles + 1))

    plt.subplot(223)
    plt.bar(range(len(rec_counts)), rec_counts)
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.xlabel('Number of recombination events')
    plt.ylabel('Molecules')
    plt.xlim((-1, num_cycles + 1))

    plt.subplot(224)
    import matplotlib.colors as mplc
    cmap = plt.get_cmap('viridis')
    cmap.set_bad(color='grey')
    plt.hist2d(cops, muts, bins=(num_cycles + 1, num_cycles + 1),
               range=[[-0.5, num_cycles + 0.5], [-0.5, num_cycles + 0.5]], norm=mplc.LogNorm(), cmap=cmap)
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))
    plt.xlabel('Number of times copied')
    plt.ylabel('Number of mutations')
    plt.gca().set_aspect('equal', 'box')

    plt.tight_layout()
    plt.savefig(f'{fp_out}', dpi=400)
