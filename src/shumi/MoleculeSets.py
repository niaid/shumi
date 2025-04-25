import logging
import random
import sys
from collections import Counter

import numpy as np
from scipy.stats import entropy, binom
from Bio.SeqIO.FastaIO import SimpleFastaParser
from . import seqlib, poplib, filelib
from .seqlib import Seq

logger = logging.getLogger(__name__)
console = logging.StreamHandler(sys.stdout)
console.setLevel(logging.INFO)
logger.addHandler(console)


class MoleculeSet:

    def __init__(self, name):
        self.name = name
        self.baseseq = None
        self.sequences = None
        self.sequences_int = None
        self.consensus = None
        self.size = None
        self.max_variant_freq = None
        self.haplotypes = None
        self.props = None
        self.counts = None
        self.entropy = None

    def write_sequences(self, fp):
        with open(fp, 'w') as f:
            for seq in self.sequences:
                f.write(f'>{seq.name}\n{seq.seq}\n')

    def read_sequences(self, fp):
        with open(fp, 'r') as f:
            self.sequences = [seqlib.Seq(name=record[0], seq=record[1]) for record in SimpleFastaParser(f)]

    def determine_consensus(self):
        if self.sequences is None:
            raise ValueError('No sequences provided before consensus calculation.')
        else:
            self.consensus = seqlib.determine_consensus(self.sequences)

    def get_entropy(self, dist=None):
        if dist is not None:
            return entropy(dist)
        if (self.props is None) and self.counts is None:
            raise ValueError('Population not initialized before entropy calculation.')
        else:
            if self.props is not None:
                self.entropy = entropy(self.props)
            else:
                self.entropy = entropy(self.counts)
        return self.entropy

    def asdict(self):
        return {seq.name: seq.seq for seq in self.sequences}


class RNASet(MoleculeSet):

    def __init__(self, name):
        super().__init__(name)
        self.mut_map = None
        logger.info(f'\tInitializing RNASet: {name}')

    def sample(self, n=1):
        return random.choices(self.sequences, weights=self.props, k=n)

    def populate_dirichlet(self, fp_base, nh, alpha):
        logger.info(f'\tUsing Dirichlet population with alpha = {alpha}, nhaps = {nh}')
        logger.info(f'\tUsing base sequence: {fp_base}')
        self.baseseq = filelib.read_fasta_seq(fp_base)
        alphas = np.array([alpha] * nh)
        self.props, self.mut_map = poplib.dirichlet_pop(self.baseseq, alphas)
        sequences = [poplib.get_cdna_seq(self.baseseq, self.mut_map, i) for i in range(nh)]
        self.sequences = []
        for i in self.mut_map:
            mut_code = self.construct_mut_code(self.mut_map[i][1])
            self.sequences.append(seqlib.Seq(name=f'RNA{i}_RNA{self.mut_map[i][0]}+{mut_code}',
                                             seq=sequences[i].seq))

    def populate_clonal(self, fp_base):
        logger.info(f'\tUsing clonal population with base sequence: {fp_base}')
        self.baseseq = filelib.read_fasta_seq(fp_base)
        self.sequences = [self.baseseq]
        self.props = [1.0]

    def populate(self, fp_pop):
        logger.info(f'\tUsing population from file: {fp_pop}')
        input_sequences = filelib.read_fasta_seqs(fp_pop)

        if not all([len(seq.seq) == len(input_sequences[0].seq) for seq in input_sequences]):
            logger.debug("All sequences provided for RNA population must be aligned to a consistent "
                         "alignment length.")
            raise ValueError("All sequences provided for RNA population must be aligned to a consistent "
                             "alignment length.")

        counts = Counter([seq.seq for seq in input_sequences])
        logger.info(f'\tNumber of sequences: {len(counts)}')

        self.sequences = [seqlib.Seq(name=f'RNA{i+1}', seq=seq) for i, seq in
                          enumerate(sorted(counts.keys(), key=lambda a: counts[a]))]
        self.props = np.array([float(counts[seq.seq]) for seq in self.sequences])
        self.props /= self.props.sum()

    def write_mut_map(self, fp):
        with open(fp, 'w') as f:
            f.write(','.join(['index', 'parent', 'mutations', 'frequency']))
            f.write('\n')
            for rna_id in self.mut_map:
                mut_code = self.construct_mut_code(self.mut_map[rna_id][1])
                f.write(','.join([str(rna_id),
                                  str(self.mut_map[rna_id][0]),
                                  mut_code,
                                  str(self.mut_map[rna_id][2])]))
                f.write('\n')

    def construct_mut_code(self, mut):
        if mut:
            mut_ind = mut[0][0]
            mut_code = f'{self.baseseq.seq[mut_ind]}{mut_ind}{mut[0][1]}'
        else:
            mut_code = 'base'
        return mut_code

    def summarize(self):
        msg = 'Summarizing RNASet:\n'
        msg += f'\tNumber of unique sequences: {len(self.sequences)}\n'
        msg += f'\tPopulation entropy: {self.get_entropy()}'
        # msg += f'Total mutation rate: {self.mu}\n'
        return msg


class FirstCopyCDNASet(MoleculeSet):

    def __init__(self, name, rnaset):
        super().__init__(name)

        self.rnaset = rnaset
        self.sampled_rna_sequences = None
        self.rna_haplotype_counts = None
        self.rna_haplotypes = None
        self.unique_cdna = None
        self.umis = None
        self.umi_count = None
        self.umi_collisions = None
        logger.info(f'\tInitializing FirstCopyCDNASet: {name}')

    def reverse_transcription(self, n, fp, rtp, fp_umi, p_err, p_recomb):

        logger.info(f'\t\tRT schema: {rtp}')
        logger.info(f'\t\tRT error rate: {p_err}')
        logger.info(f'\t\tRT recombination rate: {p_recomb}')
        if fp:
            logger.info(f'\t\tForward primer: {fp}')
        else:
            logger.debug(f'\t\tNo forward primer provided (fp = {fp})')

        umis = None
        if fp_umi is not None:
            umis = []
            logger.info(f'\tLoading reference UMIs from {fp_umi}')
            with open(fp_umi, 'r') as f:
                for line in f:
                    umi = line.strip()
                    if umi:
                        umis.append(umi)
            umis = list(set(umis))
            logger.info(f'\t\tUMI whitelist: {fp_umi} (space = {len(umis)} UMIs)')
        else:
            logger.debug(f'\t\tUsing random UMIs (umis = {fp_umi})')

        self.sampled_rna_sequences = self.rnaset.sample(n)
        self.reverse_trascribe(rtp, umis, fp=fp, p_err=p_err, p_recomb=p_recomb)

    def reverse_trascribe(self, rtp, umis, p_err, p_recomb, fp=None):

        if fp:
            fp = seqlib.sample_iupac_exclude_n(fp)
        else:
            fp = ''

        self.sequences = []
        for i, seq in enumerate(self.sampled_rna_sequences):

            recombinant = False
            # if binom.rvs(n=1, p=1 - binom.pmf(0, n=len(seq.seq), p=p_recomb)):
            if binom.rvs(n=1, p=p_recomb):

                recombinant = True
                recom_pair = self.rnaset.sample(n=1)[0]

                recomb_bp_possible = [i for i, (s1, s2) in enumerate(zip(seq.seq, recom_pair.seq)) if
                                      '-' not in (s1, s2)]

                if not recomb_bp_possible:
                    pass
                else:
                    bp_sampled = random.choice(recomb_bp_possible)
                    recombined_seq = seq.seq[:bp_sampled] + recom_pair.seq[bp_sampled:]
                    seq = Seq(name=f'recombinant+{seq.name}+{bp_sampled}+{recom_pair.name}',
                              seq=recombined_seq)

            primer, umi = seqlib.sample_primer(rtp, umis)
            core_cdna = ''.join(seqlib.int_to_nucleotide(bi) for bi in seqlib.rt(seq.int, p_err=p_err))
            cdna = seqlib.Seq(name=f'cDNA{i}_{seq.name}_{umi}',
                              seq=f'{fp + core_cdna + primer}',
                              parent=seq,
                              core=core_cdna,
                              fp=fp,
                              umi=umi,
                              rtp=primer,
                              recombinant=recombinant)
            self.sequences.append(cdna)

    def characterize(self):
        unique_rna, rna_counts = np.unique([seq.parent.seq for seq in self.sequences],
                                           return_counts=True)
        unique_cdna, cdna_counts = np.unique([seq.core for seq in self.sequences],
                                             return_counts=True)
        unique_umi, umi_counts = np.unique([seq.umi for seq in self.sequences],
                                           return_counts=True)

        self.umis = unique_umi
        self.umi_count = len(self.umis)
        self.umi_collisions = sum(umi_counts[umi_counts > 1])
        self.rna_haplotypes = unique_rna
        self.rna_haplotype_counts = rna_counts
        self.unique_cdna = unique_cdna
        self.counts = cdna_counts

    def summarize(self):
        self.characterize()
        msg = 'Summarizing FirstCopyCDNASet:\n'
        msg += f'\tNumber of sampled RNA molecules: {len(self.sequences)}\n'
        msg += f'\tNumber of unique sampled RNA sequences: {len(self.rna_haplotypes)}\n'
        msg += f'\tNumber of unique UMI used: {self.umi_count}\n'
        msg += f'\tNumber of molecules with UMI collision: {self.umi_collisions}\n'
        msg += f'\tNumber of uniquely-tagged molecules: {len(self.sequences) - self.umi_collisions}\n'
        msg += f'\tNumber of recombinant cDNA: {sum(1 if seq.recombinant else 0 for seq in self.sequences)}\n'
        msg += f'\tSampled RNA population entropy: {self.get_entropy(self.rna_haplotype_counts)}\n'
        msg += f'\tNumber of unique first-copy cDNA sequences: {len(self.unique_cdna)}\n'
        msg += f'\tFirst-copy cDNA population entropy: {self.get_entropy()}'
        # msg += f'Total mutation rate: {self.mu}\n'
        return msg


class SecondStrandSet(MoleculeSet):
    def __init__(self, name, fccdnaset, p_recomb=1e-6):
        super().__init__(name)

        self.fccdnaset = fccdnaset
        self.dumis = None
        self.dumi_count = None
        self.dumi_collisions = None
        self.p_recomb = p_recomb
        self.umis = None
        self.umi_count = None
        self.umi_collisions = None
        self.unique_cdna = None
        self.counts = None

        logger.info(f'\tInitializing SecondStrandSet: {name}')

    def second_strand_synthesis(self, fp, p_err, p_recomb, fp_umi=None):

        logger.info(f'\t\tSS-FP schema: {fp}')
        logger.info(f'\t\tSS error rate: {p_err}')
        logger.info(f'\t\tSS recombination rate: {p_recomb}')

        if not fp:
            fp = ''

        umis = None
        if fp_umi is not None:
            umis = []
            logger.info(f'\tLoading reference UMIs from {fp_umi}')
            with open(fp_umi, 'r') as f:
                for line in f:
                    umi = line.strip()
                    if umi:
                        umis.append(umi)
            logger.info(f'\t\tUMI whitelist: {fp_umi} (space = {len(umis)} UMIs)')
        else:
            logger.debug(f'\t\tUsing random UMIs (sss-umis = {fp_umi})')

        self.sequences = []
        for i, seq in enumerate(self.fccdnaset.sequences):

            recombinant = seq.recombinant

            new_seq = seq.core
            seq_name = f'ssDNA{i}_{seq.name}'

            fprimer, umi = seqlib.sample_primer(fp, umis)

            # if binom.rvs(n=1, p=1 - binom.pmf(0, n=len(seq.core), p=p_recomb)):
            if binom.rvs(n=1, p=p_recomb):
                recombinant = True
                recom_pair = random.choice(self.fccdnaset.sequences)

                recomb_bp_possible = [i for i, (s1, s2) in enumerate(zip(seq.core, recom_pair.core)) if
                                      '-' not in (s1, s2)]

                if not recomb_bp_possible:
                    pass
                else:
                    bp_sampled = random.choice(recomb_bp_possible)
                    new_seq = seq.core[:bp_sampled] + recom_pair.core[bp_sampled:]
                    seq_name = f'drecombinant+{seq.name}+{bp_sampled}+{recom_pair.name}'

            core_cdna = ''.join(seqlib.int_to_nucleotide(bi) for bi in seqlib.rt(Seq(name='', seq=new_seq).int,
                                                                                 p_err=p_err))
            ss_seq = seqlib.Seq(name=f'{seq_name}/{umi}',
                                seq=f'{fprimer + core_cdna + seq.rtp}',
                                parent=seq,
                                core=core_cdna,
                                umi=f'{seq.umi}/{umi}',
                                recombinant=recombinant)

            self.sequences.append(ss_seq)

    def characterize(self):
        unique_cdna, cdna_counts = np.unique([seq.core for seq in self.sequences],
                                             return_counts=True)
        unique_umi, umi_counts = np.unique([seq.umi for seq in self.sequences],
                                           return_counts=True)

        self.umis = unique_umi
        self.umi_count = len(self.umis)
        self.umi_collisions = sum(umi_counts[umi_counts > 1])
        self.unique_cdna = unique_cdna
        self.counts = cdna_counts

    def summarize(self):
        self.characterize()
        msg = 'Summarizing SecondStrandSet:\n'
        msg += f'\tNumber of molecules with UMI collision: {self.umi_collisions}\n'
        msg += f'\tNumber of uniquely-tagged molecules: {len(self.sequences) - self.umi_collisions}\n'
        msg += f'\tNumber of recombinant ssDNA: {sum(1 if seq.recombinant else 0 for seq in self.sequences)}\n'
        msg += f'\tNumber of unique second-strand sequences: {len(self.unique_cdna)}\n'
        msg += f'\tSecond-strand population entropy: {self.get_entropy()}'
        return msg


class PCRDNASet(MoleculeSet):

    def __init__(self, name, sequences):
        super().__init__(name)
        self.sequences = sequences
        self.intact = sum(1 for seq in self.sequences if seq.name.endswith('_intact'))
        self.recomb = sum(1 for seq in self.sequences if seq.name.endswith('_recomb'))

    def summarize(self):
        msg = 'Summarizing PCRDNASet:\n'
        msg += f'\tNumber of sampled DNA sequences: {len(self.sequences)}\n'
        msg += f'\tNumber of intact DNA sequences: {self.intact}\n'
        msg += f'\tNumber of recombinant DNA sequences: {self.recomb}'
        return msg
