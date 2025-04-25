import errno
import logging
import os
import platform
import subprocess
from Bio.SeqIO.FastaIO import SimpleFastaParser

from . import seqlib

logger = logging.getLogger(__name__)


def read_fasta_seq(fp):
    try:
        with open(fp, 'r') as f:
            name, seq = next(SimpleFastaParser(f))
    except Exception:
        logger.error('Could not read input FASTA file.')
        raise
    return seqlib.Seq(name, seq.upper())


def read_fasta_seqs(fp):
    try:
        with open(fp, 'r') as f:
            seqs = [seqlib.Seq(name=record[0], seq=record[1].upper()) for record in SimpleFastaParser(f)]
    except Exception:
        logger.error(f'Could not read input FASTA file: {fp}')
        raise
    return seqs


def write_fasta(fp, seqs):
    with open(fp, 'w') as f:
        for seq in seqs:
            f.write(f'>{seq.name}\n{seq.seq}\n')


def is_tool(name):
    try:
        devnull = open(os.devnull)
        subprocess.Popen([name], stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == errno.ENOENT:
            return False
    return True


def find_prog(prog):
    if is_tool(prog):
        cmd = "where" if platform.system() == "Windows" else "which"
        try:
            return subprocess.check_output([cmd, prog]).decode().strip()
        except subprocess.CalledProcessError:
            return None
