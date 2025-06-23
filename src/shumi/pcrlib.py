#  Copyright (c) 2025 National Institutes of Health
#  Written by Pierce Radecki
#  This program comes with ABSOLUTELY NO WARRANTY; it is intended for
#  Research Use Only and not for use in diagnostic procedures.

from dataclasses import dataclass
from typing import List, Dict, Tuple
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from functools import partial
import numpy as np
from tqdm import tqdm
import logging
from .seqlib import Seq
from .MoleculeSets import PCRDNASet

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class RecombinationEvent:
    position: int
    donor_id: str
    cycle: int


@dataclass(frozen=True)
class MutationEvent:
    position: int
    mut_ind: int
    cycle: int


@dataclass(frozen=True)
class DNAMolecule:
    parent_id: str
    molecule_id: str
    mutation_history: Tuple[RecombinationEvent | MutationEvent, ...]
    cycle_created: int
    times_copied: int
    gap_bases: Tuple[int, ...]


def get_gaps(dna_str: str) -> Tuple[int, ...]:
    return tuple(i for i, base in enumerate(dna_str) if base == '-')


class PCRSimulator:
    def __init__(self, loaded_molecules: Dict[str, str],
                 cycles: int = 30,
                 error_rate: float = 2.5e-5,
                 recombination_rate: float = 0.0001,
                 efficiency: float = 0.95,
                 max_molecules: int = 1_000_000,
                 chunk_size: int = 20_000,
                 max_threads: int = 1,
                 seed=1):
        self.initial_cdnas = loaded_molecules
        self.error_rate = error_rate
        self.num_cycles = cycles
        self.recombination_rate = recombination_rate
        self.copying_probability = efficiency
        self.max_molecules_per_cycle = max_molecules
        self.chunk_size = chunk_size
        self.max_threads = max_threads
        self.seed = seed

        self.current_pool = {
            f'{DNA_id}_':
                DNAMolecule(parent_id=DNA_id,
                            molecule_id=f'{DNA_id}_',
                            mutation_history=tuple(),
                            cycle_created=0,
                            times_copied=0,
                            gap_bases=get_gaps(loaded_molecules[DNA_id]))
            for DNA_id in loaded_molecules.keys()
        }

        self.stats_by_cycle = defaultdict(lambda: {'total_molecules': 0,
                                                   'recombined_molecules': 0})

        self.donor_table = dict()

    @staticmethod
    def _process_molecule_chunk(chunkdata: Tuple[int, List[DNAMolecule]], cycle_num: int,
                                error_rate: float,
                                recombination_rate: float,
                                copying_probability: float,
                                available_donor_ids: List[str],
                                aln_len: int,
                                root_seed: int) -> dict[str, DNAMolecule]:

        """
        Simulate a single PCR cycle on a chunk of molecules.

        Returns a dictionary of molecules that represent the PCR product.
        """

        worker_id = cycle_num * 10000 + chunkdata[0]  # Assign each worker a unique ID for rng purposes
        worker_seed = int(worker_id) + int(root_seed)
        chunk = chunkdata[1]

        rng = np.random.RandomState(worker_seed)

        new_molecules = dict()

        # Pre-populate random variables for efficiency
        muts_lookup = rng.binomial(n=aln_len, p=error_rate, size=len(chunk))
        probs = rng.random(size=(len(chunk), 2))

        for mid, molecule in enumerate(chunk):
            # mid: Molecule ID / index
            # Molecule: DNAMolecule object of molecule

            kept_molecule = DNAMolecule(parent_id=molecule.parent_id,
                                        molecule_id=f'{molecule.molecule_id}0',
                                        mutation_history=molecule.mutation_history,
                                        cycle_created=molecule.cycle_created,
                                        times_copied=molecule.times_copied,
                                        gap_bases=molecule.gap_bases)

            # Retrieve copy and recombination samples for molecule
            pcopy, precomb = probs[mid, :]

            if pcopy < copying_probability:
                # Make copy of template
                new_history = list(kept_molecule.mutation_history)

                # First simulate recombination events, to get the template bases used
                if len(available_donor_ids) > 1 and precomb < recombination_rate:
                    donor_id = None
                    while donor_id is None:
                        did = rng.choice(available_donor_ids)
                        if did != molecule.molecule_id:
                            donor_id = did

                    recomb_bp = rng.randint(1, aln_len)
                    new_history.append(RecombinationEvent(position=recomb_bp,
                                                          donor_id=donor_id,
                                                          cycle=cycle_num))

                # Simulate mutations on top of the utilized template bases
                muts = muts_lookup[mid]
                if muts:
                    pos = rng.choice(aln_len, replace=False, size=muts)  # Positions of polymerase errors
                    for i in range(muts):
                        if pos[i].item() in molecule.gap_bases:
                            continue  # If mutation in gap, we can ignore it
                        mut_ind = rng.randint(low=1, high=4)  # Uniform random mutations
                        new_history.append(MutationEvent(position=pos[i].item(), mut_ind=mut_ind, cycle=cycle_num))
                new_molecules[kept_molecule.molecule_id] = kept_molecule
                copied_mol_id = f'{molecule.molecule_id}1'

                # Save new molecule object
                new_molecules[copied_mol_id] = DNAMolecule(parent_id=molecule.parent_id,
                                                           molecule_id=copied_mol_id,
                                                           mutation_history=tuple(new_history),
                                                           cycle_created=cycle_num,
                                                           times_copied=molecule.times_copied + 1,
                                                           gap_bases=molecule.gap_bases)
            else:
                new_molecules[kept_molecule.molecule_id] = kept_molecule

        return new_molecules

    def simulate(self) -> Dict:

        sequence_lengths = {seq_id: len(seq) for seq_id, seq in self.initial_cdnas.items()}
        aln_len = next(iter(sequence_lengths.values()))

        for cycle_num in range(1, self.num_cycles + 1):

            working_pool = list(self.current_pool.values())

            # Down-sample PCR pool if needed
            if len(working_pool) > self.max_molecules_per_cycle:
                working_pool = np.random.choice(working_pool, size=self.max_molecules_per_cycle, replace=False)

            # Make chunks of molecules to simulate cycle in parallel
            chunk_generator = chunkify(working_pool, self.chunk_size)

            # List of molecule IDs to serve as putative recombinant partners
            available_donor_ids = [mol.molecule_id for mol in working_pool]

            num_chunks = int(np.ceil(len(working_pool) / self.chunk_size))
            new_pool = dict()
            num_processes = min(self.max_threads, max(1, len(working_pool) // self.chunk_size))

            # Progress bar that launches parallel jobs
            with tqdm(total=num_chunks, leave=False, disable=True if num_chunks < 2 else False, miniters=1,
                      unit=' chunks', bar_format=f'        PCR cycle {cycle_num}' + '{l_bar}{bar:10}{r_bar}') as pbar:
                with ProcessPoolExecutor(max_workers=num_processes) as executor:
                    process_func = partial(self._process_molecule_chunk,
                                           error_rate=self.error_rate,
                                           cycle_num=cycle_num,
                                           recombination_rate=self.recombination_rate,
                                           available_donor_ids=available_donor_ids,
                                           aln_len=aln_len,
                                           copying_probability=self.copying_probability,
                                           root_seed=self.seed)

                    # Simulate PCR in chunks
                    chunk_results = executor.map(process_func, chunk_generator)
                    for result in chunk_results:
                        pbar.update()
                        new_pool.update(result)  # Save new molecules returned by each chunk

            # After cycle, track which donors were involved in recombination to track separately
            donors = [e.donor_id for mol in new_pool.values() for e in mol.mutation_history if
                      type(e) is RecombinationEvent]
            self.donor_table.update(
                {donor: self.current_pool[donor] for donor in donors if donor not in self.donor_table})

            logger.debug(f'\t\tNumber of molecules in PCR pool after cycle {cycle_num}: {len(new_pool)}')

            self.current_pool = new_pool
            self._update_stats(cycle_num)

        return self.stats_by_cycle

    def _update_stats(self, cycle_num: int) -> None:

        stats = self.stats_by_cycle[cycle_num]
        stats['total_molecules'] = len(self.current_pool)
        recombined = sum(1 for mol in self.current_pool.values() if
                         any(type(mut) is RecombinationEvent for mut in mol.mutation_history))
        stats['recombined_molecules'] = recombined
        logger.debug(f'\t\t\tNumber of PCR recombinants after cycle {cycle_num}: {recombined}'
                     f' ({recombined / len(self.current_pool) * 100:1.4f}%)')

    def reconstruct_sequence(self, molecule: DNAMolecule):
        """
        Reconstruct the nucleotide sequence from a DNA molecule using its mutation history.
        """
        if not molecule.mutation_history:
            return self.initial_cdnas[molecule.parent_id]

        sequence = self.initial_cdnas[molecule.parent_id]
        for event in molecule.mutation_history:

            # Handle recombination event
            if type(event) is RecombinationEvent:
                bp = event.position
                donor = event.donor_id
                if donor in self.initial_cdnas:
                    sequence = self.initial_cdnas[donor][:bp] + sequence[bp:]
                else:
                    sequence = self.reconstruct_sequence(self.donor_table[donor])[:bp] + sequence[bp:]

            # Handle mutations
            if type(event) is MutationEvent:
                bp = event.position
                mut = event.mut_ind
                sequence = sequence[:bp] + apply_mut(sequence[bp], mut) + sequence[bp + 1:]

        return sequence

    def sample(self, nseqs):
        """
        Sample a number of molecules from the PCR product. Returns PCRDNASet.
        """

        nseqs = int(nseqs)
        try:
            seq_names = np.random.choice(list(self.current_pool.keys()), size=nseqs, replace=False)
        except ValueError:
            print(f'Not enough sequences in pool for sampling desired number of molecules.\n'
                  f'Number of sequences in pool: {len(self.current_pool)}\n'
                  f'Number of molecules requested: {nseqs}')
            raise
        seqs = []
        for seq in seq_names:
            seq_name = seq + '_recomb' if \
                any(type(e) is RecombinationEvent for e in self.current_pool[seq].mutation_history) \
                else seq + '_intact'
            seqs.append(Seq(name=seq_name, seq=self.reconstruct_sequence(self.current_pool[seq]),
                            parent=self.current_pool[seq].parent_id))
        return PCRDNASet(name='pcr_cdna', sequences=seqs)


nuc2ind = {'A': 0,
           'C': 1,
           'G': 2,
           'T': 3,
           '-': 4}

ind2nuc = {ind: base for base, ind in nuc2ind.items()}


def apply_mut(base, mut_ind):
    """Apply a sampled mutation to a nucleotide."""
    if base == '-':
        return '-'
    bi = nuc2ind[base]
    bi += mut_ind
    bi %= 4
    return ind2nuc[bi]


def chunkify(lst, n):
    """Yield successive chunks from array."""
    for i in range(0, len(lst), n):
        yield i + 1, lst[i:i + n]
