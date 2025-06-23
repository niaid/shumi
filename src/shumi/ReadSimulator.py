#  Copyright (c) 2025 National Institutes of Health
#  Written by Pierce Radecki
#  This program comes with ABSOLUTELY NO WARRANTY; it is intended for
#  Research Use Only and not for use in diagnostic procedures.

import functools
import logging
import multiprocessing
import numpy as np
import os
import subprocess
from pathlib import Path
import pysam
from matplotlib import pyplot as plt
from tqdm import tqdm
from scipy.stats import nbinom

from . import filelib, seqlib

logger = logging.getLogger(__name__)


# noinspection PyProtectedMember
def clipped_nb(p):
    """
    Computes PMF for a clipped negative binomial distribution with support [3, 60].
    """

    pmf = np.array([nbinom._pmf(i, 2, p) for i in range(61)])
    pmf[0] = 0.0
    pmf[1] = 0.0
    pmf[2] = 0.0
    pmf[60] += nbinom._sf(61, 2, p)
    return pmf / pmf.sum()


def get_slope(mt, sp1, sp2):
    return sp1 * mt + sp2
    # return -0.00000102 * mt + 0.00004281


def get_nbin_p(rl, mt, sp1=-0.00000102, sp2=0.00004281):
    slope = get_slope(mt, sp1, sp2)
    return slope * rl


def get_dist(rl, mt, sp1, sp2):
    return clipped_nb(get_nbin_p(rl, mt, sp1, sp2))


def num_passes(cdna_length, movie_time, sp1, sp2):
    dist_sample = get_dist(cdna_length, movie_time, sp1, sp2)
    return np.random.choice(np.arange(61), 1, p=dist_sample)[0]


class ReadSimulator:
    def __init__(self, name, reads, run_dir, threads=1, output=None, pcrdna_set=None, movietime=20,
                 sp1=-0.00000102, sp2=0.00004281, rng=1, ccs=None):

        self.name = name
        self.reads = reads
        self.threads = threads
        self.output = output
        self.pcrdna_set = pcrdna_set
        self.movietime = movietime
        self.sp1 = sp1
        self.sp2 = sp2
        self.rng = rng
        self.passes_dict = dict()
        self.output_dir = Path(self.output) if self.output is not None else Path.cwd()
        self.run_dir = run_dir
        self.ccs_produced = None
        self.reads_produced = None
        self.ccs_name2pcr_name = None

        if ccs is None:
            if filelib.find_prog('ccs'):
                logger.debug('ccs command found, will form CCS reads.')
                self.ccs = True
            else:
                logger.debug('ccs command not found, will not form CCS reads.')
                self.ccs = False
        else:
            self.ccs = ccs

    def set_pcrcdna_set(self, pcrcdna_set):
        self.pcrdna_set = pcrcdna_set

    def sample_num_passes(self):
        for seq in self.pcrdna_set.sequences:
            self.passes_dict[seq.name] = num_passes(seq.true_len, self.movietime, self.sp1, self.sp2)
        return

    def split_by_num_passes(self):
        os.makedirs(self.output_dir, exist_ok=True)
        for passes in set(self.passes_dict.values()):
            seqs = [pcrcdna for pcrcdna in
                    self.pcrdna_set.sequences
                    if self.passes_dict[pcrcdna.name] == passes]

            fp_num_passes = os.path.join(self.output_dir, f'{self.name}.{passes}-passes.fasta')
            with open(fp_num_passes, 'w') as f:
                for seq in seqs:
                    f.write(f'>{seq.name}_{passes}-passes\n{seq.seq.replace("-", "")}\n')

        return

    def get_mean_passes(self):
        return sum([x * sum(self.passes_dict[i] == x for i in self.passes_dict) / len(self.passes_dict)
                    for x in set(self.passes_dict.values())])

    def summarize(self):
        msg = 'Summarizing ReadSimulator:\n'
        msg += f'\tNumber of target reads: {self.reads}\n'
        msg += f'\tNumber of molecules sequenced: {len(self.pcrdna_set.sequences)}\n'
        msg += f'\tMean number of passes: {self.get_mean_passes():.2f}\n'
        msg += f'\tNumber of CCS produced: {self.ccs_produced}\n'
        msg += f'\tNumber of CCS sampled: {self.reads_produced}'
        return msg

    def get_reads(self, pcrcdnaset=None):
        if pcrcdnaset is not None:
            self.pcrdna_set = pcrcdnaset
        else:
            if not self.pcrdna_set:
                raise ValueError("Provide a PCRDNASet before getting reads. ")

        self.sample_num_passes()
        self.split_by_num_passes()

        fig = plt.figure(figsize=(6, 5))
        ax = fig.add_subplot(111)
        ax.hist(self.passes_dict.values(), bins=np.arange(61))
        ax.set_xlabel('Number of passes')
        ax.set_ylabel('Count')
        plt.savefig(os.path.join(self.run_dir, f'{self.name}-passes.png'), format='png')

        pool = multiprocessing.Pool(processes=self.threads)

        logger.debug('\tGenerating subreads...')

        npasses_iter = sorted(set(self.passes_dict.values()), reverse=True)

        with tqdm(total=len(npasses_iter), leave=False, unit=' batches') as pbar:
            try:
                run_func = functools.partial(self.generate_subreads_for_passes_count,
                                             output_dir=self.output_dir,
                                             name=self.name,
                                             rng=self.rng)
                for _ in pool.imap_unordered(run_func, npasses_iter):
                    pbar.update()
                pool.close()
                pool.join()
            except Exception:
                pool.terminate()
                raise

        # Merge subreads from each number-of-passes count
        logger.info('\tMerging subread .bam files...')
        cmd = ' '.join(['samtools', 'merge', '-r', '-f',
                        f'{self.run_dir}/{self.name}.subreads-merged.bam',  # Output file name
                        f'{self.output_dir}/*.bam'])  # Input files
        p = subprocess.Popen(cmd, shell=True)  # shell = True necessary for globbing
        p.communicate()

        cmd = ' '.join(['samtools', 'sort', '-n', '-o', f'{self.run_dir}/{self.name}.subreads.bam',
                        f'{self.run_dir}/{self.name}.subreads-merged.bam'])  # Input files
        p = subprocess.Popen(cmd, shell=True)
        p.communicate()

        # Generate CCS, if command is available
        if self.ccs:
            self.generate_ccs()
        else:
            logger.warning('\tSkipping CCS generation as the ccs/pbccs command was not configured.')
            self.ccs_produced = 0
            self.reads_produced = 0

        # Compress number of passes directory containing intermediate files
        cmd = ['tar', '-czf', f'{self.output_dir}.tar.gz',
               '-C', str(self.output_dir), '.']
        p = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        p.communicate()

        # Remove uncompressed number of passes directory
        cmd = ['rm', '-r', str(self.output_dir)]
        p = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        p.communicate()

        return

    @staticmethod
    def generate_subreads_for_passes_count(passes, output_dir, name, rng):

        pbsim_src = Path(filelib.find_prog('pbsim')).parent
        data_model = os.path.join(pbsim_src, '../data/ERRHMM-SEQUEL.model')
        fp_num_passes = os.path.join(output_dir, f'{name}.{passes}-passes.fasta')

        with open(f'{output_dir}/{name}.{passes}-passes.err', 'w') as ferr, \
                open(f'{output_dir}/{name}.{passes}-passes.out', 'w') as fout:
            cmd = ['pbsim', '--seed', f'{rng}', '--strategy', 'templ',
                   '--method', 'errhmm',
                   '--errhmm', f'{data_model}',
                   '--template', f'{fp_num_passes}',
                   '--pass-num', f'{passes}',
                   '--difference-ratio', '22:45:33',
                   '--prefix', f'{output_dir}/{name}.{passes}-passes',
                   '--id-prefix', f'{name}.{passes}',
                   ]
            subprocess.call(cmd, stdout=fout, stderr=ferr)

        logger.info('\tCompresing subreads')
        cmd = ['samtools', 'view', '-S', '-b',
               '-o', f'{output_dir}/{name}.{passes}-passes.bam',
               f'{output_dir}/{name}.{passes}-passes.sam']
        p = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        p.communicate()

        cmd = ['rm', f'{output_dir}/{name}.{passes}-passes.sam',
               f'{output_dir}/{name}.{passes}-passes.maf']
        p = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        p.communicate()

    def generate_ccs(self):
        logger.info('\tForming ccs reads...')
        with open(f'{self.run_dir}/{self.name}.ccs.out', 'w') as f_out:
            cmd = ['ccs', '--num-threads', f'{self.threads}', '--log-level', 'INFO',
                   f'{self.run_dir}/{self.name}.subreads.bam', f'{self.run_dir}/{self.name}.ccs_prelim.bam']
            p = subprocess.Popen(cmd, stdout=f_out, stderr=f_out)
        p.communicate()

        logger.info('\tDownsampling reads...')
        bam_ccs = pysam.AlignmentFile(f'{self.run_dir}/{self.name}.ccs_prelim.bam', mode='rb', check_sq=False)

        ccs_read_ids = []
        for read in bam_ccs:
            ccs_read_ids.append(read.query_name)

        ccs_read_ids.sort()

        self.ccs_produced = len(ccs_read_ids)

        if self.reads > len(ccs_read_ids):
            error_msg = f"""
Fatal error: insufficient number of CCS reads generated. Increase --overhead to ensure
enough CCS are generated for the run.\n\n
Requested reads:\t{self.reads}\n
Molecules sequenced:\t{len(self.pcrdna_set.sequences)}\n
Produced reads:\t{self.ccs_produced}\n\nRun failed.\n
            """
            logger.error(error_msg)
            raise ValueError(error_msg.replace('\n', ' '))

        sampled_reads = np.random.choice(ccs_read_ids, size=self.reads, replace=False)

        with open(f'{self.run_dir}/{self.name}.ccs_sampled.txt', 'w') as f:
            f.write('\n'.join(sampled_reads))

        with open(f'{self.run_dir}/{self.name}.ccs.bam', 'w') as f:
            cmd = ['samtools', 'view', '-b',
                   '-N', f'{self.run_dir}/{self.name}.ccs_sampled.txt',
                   f'{self.run_dir}/{self.name}.ccs_prelim.bam']
            subprocess.call(cmd, stdout=f)

        ccs_np_dict = dict()
        for read in sampled_reads:
            passes = int(read.split('/')[0].split('.')[-1])
            # ind = int(read.split('/')[1])
            if passes in ccs_np_dict:
                ccs_np_dict[passes].append(read)
            else:
                ccs_np_dict[passes] = [read]

        self.ccs_name2pcr_name = dict()

        for passes in ccs_np_dict:
            passes_read_names = seqlib.get_read_names(
                os.path.join(self.output_dir, f'{self.name}.{passes}-passes.fasta'))
            for read in ccs_np_dict[passes]:
                passes_ind = int(read.split('/')[1])
                self.ccs_name2pcr_name[read] = '_'.join(passes_read_names[passes_ind - 1].split('_')[:-1])

        with open(f'{self.run_dir}/{self.name}.ccs2cdna.txt', 'w') as f:
            for read in sorted(sampled_reads, key=lambda a: (int(a.split('/')[0].split('.')[-1]),
                                                             int(a.split('/')[1]))):
                f.write(f'{read}\t{self.ccs_name2pcr_name[read]}\n')
        self.reads_produced = len(sampled_reads)

        logger.debug('\tWriting FASTQ output...')
        cmd = ['samtools', 'bam2fq', '-T', 'np,rq',
               '-0', f'{self.run_dir}/{self.name}.ccs.fastq.gz',
               f'{self.run_dir}/{self.name}.ccs.bam']
        p = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        p.communicate()

        logger.debug('\tWriting FASTA output...')
        cmd = ['seqkit', 'fq2fa', '-o', f'{self.run_dir}/{self.name}.ccs.fasta',
               f'{self.run_dir}/{self.name}.ccs.fastq.gz']
        p = subprocess.Popen(cmd)
        p.communicate()

        logger.info('\tdone')
