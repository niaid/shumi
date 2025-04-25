import logging
import os
import random
import sys
from collections import Counter
from random import seed

import numpy as np

from . import arglib
from . import filelib
from . import graphlib
from . import seqlib
from . import version
from .MoleculeSets import RNASet, FirstCopyCDNASet, SecondStrandSet
from .ReadSimulator import ReadSimulator
from .pcrlib import PCRSimulator


def main(test_command=None):
    if test_command is None:
        args = arglib.parse_args(sys.argv[1:])
    else:
        args = arglib.parse_args(test_command.split())

    if args.rng:
        rng = int(args.rng)
        seed(rng)
        np.random.seed(rng)
    else:
        rng = random.randint(100000000, 1000000000 - 1)
        seed(rng)
        np.random.seed(rng)

    os.makedirs(args.output, exist_ok=True)

    # Configure logger
    logging.basicConfig(filename=os.path.join(args.output, f'{args.name}-{args.hex_id}.log'),
                        format='%(asctime)s %(levelname)s %(name)s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.DEBUG)

    # Split logging to file (DEBUG) and stdout (INFO)
    logger = logging.getLogger(__name__)
    logging.getLogger('matplotlib').setLevel(logging.WARN)
    # logger.setLevel(logging.DEBUG)

    # Configure stdout messages
    console = logging.StreamHandler(sys.stdout)
    console.setLevel(logging.INFO)
    logger.addHandler(console)

    logger.info(graphlib.logo)
    logger.info(f'shumi version {version.__version__}')
    logger.info(f'Args: {arglib.print_selected_args(args, args.__dict__.keys())}')

    logger.info(f'Using RNG seed: {rng}')
    logger.debug('test')

    run_name = args.name

    logger.info('Initializing RNASet...')
    rna_set = RNASet(f'{run_name}_rna')

    if args.random:
        args.input_seq = os.path.join(args.output, f'{args.name}-base-{args.hex_id}.fasta')
        rseq = seqlib.random_seq(args.random)
        base_seq = seqlib.Seq(name=f'random-base-seq-{args.hex_id}', seq=rseq)
        filelib.write_fasta(args.input_seq, [base_seq])

    if args.input_seq:
        if args.dirichlet:
            rna_set.populate_dirichlet(args.input_seq, args.nhaps, args.dirichlet)
            rna_set.write_mut_map(os.path.join(args.output, f'{args.name}-mut-map.txt'))
        else:
            rna_set.populate_clonal(args.input_seq)

    elif args.input_pop:
        rna_set.populate(args.input_pop)

    rna_set.write_sequences(os.path.join(args.output, f'{args.name}-rna.fasta'))
    logger.info('\tdone')
    logger.info(rna_set.summarize())

    logger.info('Initializing FirstCopyCDNASet...')
    fccdna_set = FirstCopyCDNASet(name=f'{run_name}_cdna', rnaset=rna_set)
    fccdna_set.reverse_transcription(n=args.ncdna, fp=args.fp, rtp=args.rtp, fp_umi=args.umis,
                                     p_err=args.error_rate_rt, p_recomb=args.recomb_rate_rt)
    fccdna_set.write_sequences(os.path.join(args.output, f'{args.name}-fccdna.fasta'))
    logger.info('\tdone')
    logger.info(fccdna_set.summarize())

    if args.sss:
        logger.info(f'Initializing SecondStrandSet...')
        ss_set = SecondStrandSet(name=f'{run_name}_ssdna', fccdnaset=fccdna_set)
        ss_set.second_strand_synthesis(args.sss_fp, fp_umi=args.sss_umis,
                                       p_err=args.error_rate_sss, p_recomb=args.recomb_rate_sss)
        ss_set.write_sequences(os.path.join(args.output, f'{args.name}-ssdna.fasta'))
        logger.info('\tdone')
        logger.info(ss_set.summarize())
        pcr_input = ss_set
    else:
        pcr_input = fccdna_set

    logger.info('Simulating PCR reaction...')

    pcr_simulator = PCRSimulator(cycles=args.cycles,
                                 efficiency=args.pcr_efficiency,
                                 error_rate=args.error_rate_pcr,
                                 recombination_rate=args.recomb_rate_pcr,
                                 max_molecules=args.max_molecules_tracked,
                                 loaded_molecules=pcr_input.asdict(),
                                 max_threads=args.threads,
                                 seed=rng)

    pcr_simulator.simulate()
    pcrdna_set = pcr_simulator.sample(np.ceil(args.reads * args.overhead))

    logger.info('\tdone')
    logger.info(pcrdna_set.summarize())

    logger.info('Writing sampled PCR-derived DNA sequences...')
    pcrdna_set.write_sequences(os.path.join(args.output, f'{args.name}-pcrdna.fasta'))
    logger.info('\tdone')

    graphlib.pcr_summary_graph(pcr_simulator, os.path.join(args.output, f'{args.name}-pcr-summary.png'))

    logger.info('Simulating sequencing of PCR product...')
    read_simulator = ReadSimulator(reads=args.reads,
                                   pcrdna_set=pcrdna_set,
                                   run_dir=args.output,
                                   output=os.path.join(args.output, f'{run_name}-{args.hex_id}-passes/'),
                                   name=run_name,
                                   movietime=args.movie_time,
                                   sp1=args.sp1,
                                   sp2=args.sp2,
                                   rng=rng,
                                   threads=args.threads)
    read_simulator.get_reads()

    logger.info('\tdone')
    logger.info(read_simulator.summarize())

    logger.info('Building summary graphs...')

    # Count molecules from each parent cDNA in PCR output
    pcr_cdna = Counter([mol.parent_id for mol in pcr_simulator.current_pool.values()])
    pcr_counts = [pcr_cdna[seq.name] for seq in pcr_input.sequences]

    # Count molecules from each parent cDNA sampled for sequencing
    sampled_cdna = Counter([seq.parent for seq in pcrdna_set.sequences])
    sampled = [sampled_cdna[seq.name] for seq in pcr_input.sequences]

    # Count CCS from each parent cDNA
    if read_simulator.ccs_name2pcr_name:
        ccs_cdna = Counter(
            [seq.parent for seq in pcrdna_set.sequences if seq.name in read_simulator.ccs_name2pcr_name.values()])
        sampled_ccs = [ccs_cdna[seq.name] for seq in pcr_input.sequences]
    else:
        sampled_ccs = None

    graphlib.summary_graph(pcr_counts, sampled, sampled_ccs,
                           fp_out=os.path.join(args.output, f'{args.name}-samples.png'))
    logger.info('\tdone')

    logger.info(f'shumi run {args.name}-{args.hex_id} complete.\n')


if __name__ == "__main__":
    main()
