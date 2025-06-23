#  Copyright (c) 2025 National Institutes of Health
#  Written by Pierce Radecki
#  This program comes with ABSOLUTELY NO WARRANTY; it is intended for
#  Research Use Only and not for use in diagnostic procedures.

import numpy as np
import random
from Bio.SeqIO.FastaIO import SimpleFastaParser

mut_opts = {
    'A': ['C', 'G', 'T'],
    'C': ['A', 'G', 'T'],
    'G': ['A', 'C', 'T'],
    'T': ['A', 'C', 'G'],
    '-': ['-']
}

iupac_ambiguous_to_unambiguous = {
    'A': ['A'],
    'C': ['C'],
    'G': ['G'],
    'T': ['T'],
    'U': ['U'],
    'R': ['A', 'G'],
    'Y': ['C', 'T'],
    'S': ['G', 'C'],
    'W': ['A', 'T'],
    'K': ['G', 'T'],
    'M': ['A', 'C'],
    'B': ['C', 'G', 'T'],
    'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'],
    'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T']
}

iupac_ambiguous_to_unambiguous_exclude_n = {
    'N': ['N'],
    'A': ['A'],
    'C': ['C'],
    'G': ['G'],
    'T': ['T'],
    'U': ['U'],
    'R': ['A', 'G'],
    'Y': ['C', 'T'],
    'S': ['G', 'C'],
    'W': ['A', 'T'],
    'K': ['G', 'T'],
    'M': ['A', 'C'],
    'B': ['C', 'G', 'T'],
    'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'],
    'V': ['A', 'C', 'G']
}

mut_rates = {
    'A': [1 - 1e-4, 1e-4 / 3, 1e-4 / 3, 1e-4 / 3],
    'C': [1e-4 / 3, 1 - 1e-4, 1e-4 / 3, 1e-4 / 3],
    'T': [1e-4 / 3, 1e-4 / 3, 1 - 1e-4, 1e-4 / 3],
    'G': [1e-4 / 3, 1e-4 / 3, 1e-4 / 3, 1 - 1e-4]
}


def parse_fasta_seq(fp_fasta):
    try:
        with open(fp_fasta, 'r') as f:
            return next(SimpleFastaParser(f))[1].upper()
    except Exception as e:
        print(f'Could not parse fasta file: {fp_fasta}')
        raise e


def get_read_names(fp_fasta):
    with open(fp_fasta, 'r') as f:
        seqs = [record[0] for record in SimpleFastaParser(f)]
    return seqs


def sample_primer(primer, umis=None):

    primer = sample_iupac_exclude_n(primer)

    pout = ''
    umi_tag = None

    if primer.count('N') == 0:
        raise ValueError(f'No UMI region found in RT primer: {primer}')

    if umis is None:
        for p in primer:
            if p == 'N':
                pout += random.choice(['A', 'C', 'T', 'G'])
            else:
                pout += p
    else:
        segments = [seg for seg in primer.split('N') if seg]
        umis_primer = []
        pout = segments[0]
        for seg in segments[1:]:
            umi = random.choice(umis)
            umis_primer.append(umi)
            pout += umi + seg
        umi_tag = '-'.join(umis_primer)

    if umi_tag:
        pass
    else:
        umi_tag = ''.join([pout[i] for i in range(len(pout)) if primer[i] == 'N'])

    return pout, umi_tag


def sample_iupac_exclude_n(primer):
    return ''.join([random.choice(iupac_ambiguous_to_unambiguous_exclude_n[b]) for b in primer])


def rt(seq, p_err=1e-4):
    mut_pos = np.random.binomial(1, p_err, size=len(seq))
    mut_pos[mut_pos > 0] *= np.random.randint(low=1, high=3, size=len(mut_pos[mut_pos > 0]))
    rt_seq = (seq + mut_pos) % 4
    rt_seq[seq == 4] = 4
    return rt_seq


mappingbase = {'A': 0, 'C': 1, 'G': 2, 'T': 3, '-': 4}
mappingint = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: '-'}


def nucleotide_to_int(nucleotide):
    return mappingbase[nucleotide]


def int_to_nucleotide(bi):
    return mappingint[bi]


mut_base_options = {
    'A': ['T', 'G', 'C'],
    'T': ['A', 'G', 'C'],
    'G': ['A', 'T', 'C'],
    'C': ['A', 'T', 'G'],
    '-': ['-']
}

mut_int_options = {nucleotide_to_int(i): [nucleotide_to_int(j) for j in mut_base_options[i]] for i in mut_base_options}


def determine_consensus(seqs):
    consensus = ''
    for i in range(len(seqs[0].seq)):
        bases = [seq.seq[i] for seq in seqs]
        counts = {base: bases.count(base) for base in bases}
        consensus += max(counts, key=counts.get)
    consensus_seq = Seq('Consensus', consensus)
    return consensus_seq


def random_seq(nbases):
    return ''.join(random.choices(['A', 'C', 'G', 'T'], k=nbases))


class Seq:
    def __init__(self, name, seq, parent=None, core=None, fp=None, umi=None, rtp=None, recombinant=False):
        self.name = str(name)
        self.seq = str(seq)
        self.int = np.array([nucleotide_to_int(base) for base in self.seq], dtype=np.int8)
        self.len = len(seq)
        self.true_len = len(seq) - seq.count('-')
        self.parent = parent
        self.core = core
        self.fp = fp
        self.umi = umi
        self.rtp = rtp
        self.recombinant = recombinant
        self.pcr_counts = None
        self.pcr_samples = None

    def __len__(self):
        return self.len
