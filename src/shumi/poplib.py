#  Copyright (c) 2025 National Institutes of Health
#  Written by Pierce Radecki
#  This program comes with ABSOLUTELY NO WARRANTY; it is intended for
#  Research Use Only and not for use in diagnostic procedures.

import logging
import numpy as np
from scipy.stats import dirichlet

from . import seqlib

logger = logging.getLogger(__name__)


def dirichlet_pop(base_seq, alphas):
    logger.info(f'\tGenerating Dirichlet RNA population, alpha={alphas[0]}')

    alphas = np.array(alphas)

    dist = dirichlet(alpha=alphas)
    dist_sample = np.array([0.0] * len(alphas))

    while dist_sample.sum() != 1.0:
        dist_sample = np.array(dist.rvs())[0]
        dist_sample[dist_sample <= 0] = 1e-280
        dist_sample[np.isnan(dist_sample)] = 1e-280
        dist_sample.sort()
        dist_sample = dist_sample[::-1]
        dist_sample /= dist_sample.sum()

    logger.info(f'\tDirichlet distribution sampled')
    mut_map = {}

    sampled_muts = np.arange(len(base_seq))
    np.random.shuffle(sampled_muts)

    for i in range(len(alphas)):
        if i == 0:
            mut_map[0] = (0, [], dist_sample[i])
        else:
            pseq = np.random.choice(i, p=dist_sample[:i] / dist_sample[:i].sum())
            mut, sampled_muts = sampled_muts[-1], sampled_muts[:-1]
            mut_base = np.random.choice(seqlib.mut_opts[base_seq.seq[mut]])
            mut_map[i] = (pseq, [(mut, mut_base)], dist_sample[i])

    logger.info(f'\tMutation tree sampled')

    return dist_sample, mut_map


def get_cdna_seq(base_seq, mm, rna_id):
    if not mm[rna_id][1]:
        return base_seq
    else:
        pseq = get_cdna_seq(base_seq, mm, mm[rna_id][0]).seq
        mut_ind = mm[rna_id][1][0][0]
        new_base = mm[rna_id][1][0][1]
        return seqlib.Seq(name='', seq=pseq[:mut_ind] + new_base + pseq[mut_ind + 1:])

