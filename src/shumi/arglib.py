import argparse
import datetime
import random
from pathlib import Path
import yaml
from . import version


def parse_args(args):
    parser = argparse.ArgumentParser(description="A command line tool to simulate HT-SGS data.",
                                     usage="%(prog)s [--input-seq path | --input-pop path | --random L] "
                                           "<--ncdna n> <--rtp ACGTN...> <--reads m> [-h] [options]")

    misc_parser = parser.add_argument_group("Miscellaneous Options")
    misc_parser.add_argument("--config", metavar='path', default=None, type=str, required=False,
                             help="Use arguments in config file. Command line arguments override config.")
    misc_parser.add_argument("--output", metavar="path", type=str, default=None, help="Output directory")
    misc_parser.add_argument("--name", metavar="path", type=str, default=None,
                             help="Run name for output files")
    misc_parser.add_argument("--rng", metavar='digits', default=None, type=int, required=False,
                             help="Use input seed for reproducibility")
    misc_parser.add_argument("--threads", metavar="t", type=int, default=1,
                             help="Number of processes to use")
    misc_parser.add_argument("-v", "--version", action='version', help="Print version number and exit",
                             version='shumi {version}'.format(version=version.__version__))

    rnapop_parser = parser.add_argument_group("RNA Population Options")
    rnapop_parser.add_argument("--input-seq", metavar="path", type=str, required=False,
                               help="Input file path to base sequence for simulation")
    rnapop_parser.add_argument("--input-pop", metavar="path", type=str, required=False,
                               help="Input file path to reference population for simulation")
    rnapop_parser.add_argument("--dirichlet", metavar='alpha', type=float, required=False,
                               help="Use dirichlet model of haplotype frequencies in RNA population")
    rnapop_parser.add_argument("--nhaps", metavar='h', type=int, required=False,
                               help="Number of haplotypes considered in dirichlet population model")
    rnapop_parser.add_argument("--random", metavar='L', type=int, required=False,
                               help="Use random sequence of length L to simulate RNA population")
    rnapop_parser.add_argument("--clonal", action='store_true', required=False,
                               help="Use clonal RNA population")

    cdna_parser = parser.add_argument_group("First-Copy cDNA Options")
    cdna_parser.add_argument("--ncdna", metavar='n', type=int, default=None, required=False,
                             help="** Number of first-copy cDNA to be generated")
    cdna_parser.add_argument("--rtp", metavar='[ACGTN...]', type=str, default=None, required=False,
                             help="** RT primer sequence. Runs of N will replaced at random.")
    cdna_parser.add_argument("--fp", metavar='[ACGT...]', type=str, default=None,
                             help="Optional sequence to append to 5' end of template sequence during RT")
    cdna_parser.add_argument("--umis", metavar='path', type=str, default=None,
                             help="File of reference UMIs to use when generating cDNA")
    cdna_parser.add_argument("--error-rate-rt", metavar='p', type=float, default=1e-4,
                             help="Uniform per-base error rate for reverse transcription (default: 0.0001)")
    cdna_parser.add_argument("--recomb-rate-rt", metavar='xr', type=float, default=0.005,
                             help="Probability of recombination when copying a molecule during RT (default: 0.005)")

    ss_parser = parser.add_argument_group("Second Strand Synthesis Options")
    ss_parser.add_argument("--sss", action='store_true', default=False, required=False,
                           help="Perform second strand synthesis with UMI")
    ss_parser.add_argument("--sss-fp", metavar='[ACGTN...]', type=str, default=None, required=False,
                           help="Second strand synthesis primer sequence. Runs of N will replaced at random.")
    ss_parser.add_argument("--sss-umis", metavar='path', type=str, default=None,
                           help="File of reference UMIs to use when peforming second strand synthesis")
    ss_parser.add_argument("--error-rate-sss", metavar='p', type=float, default=None,
                           help="Uniform per-base error rate for second strand synthesis (default: same as PCR)")
    ss_parser.add_argument("--recomb-rate-sss", metavar='xr', type=float, default=None,
                           help="Probability of recombination when copying a molecule during second strand synthesis "
                                "(default: same as PCR)")

    pcr_parser = parser.add_argument_group("PCR Options")
    pcr_parser.add_argument("--cycles", metavar='k', type=int, default=30, help="Number of PCR cycles "
                                                                                "(default: 30)")
    pcr_parser.add_argument("--pcr-efficiency", metavar='q', type=float, default=0.95,
                            help="Efficiency of PCR simulation; i.e., probability a molcule is copied during a cycle "
                                 "(default: 0.95)")
    pcr_parser.add_argument("--error-rate-pcr", metavar='r', type=float, default=2.6e-5,
                            help="Uniform per-base error rate for PCR copies (default: 2.5e-5)")
    pcr_parser.add_argument("--recomb-rate-pcr", metavar='s', type=float, default=0.0001,
                            help="Probability of recombination when copying a molecule during PCR (default: 0.0001)")
    pcr_parser.add_argument("--max-molecules-tracked", metavar='b', type=int, default=1000000,
                            help="Maximum number of molecules to track during PCR cycles (default: 1000000)")

    pbseq_parser = parser.add_argument_group("Sequencing Options")
    pbseq_parser.add_argument("--reads", metavar='m', type=int, default=None, required=False,
                              help="** Desired number of CCS reads")
    pbseq_parser.add_argument("--overhead", metavar='u', type=float, default=1.5,
                              help="Ratio for number PCR DNA molecules to generate subreads for (default: 1.5)")
    pbseq_parser.add_argument("--movie-time", metavar='mt', type=float, default=20.0, required=False,
                              help="Movie time of sequencing in hours, must be between 10 and 40 (default: 20)")
    pbseq_parser.add_argument("--sp1", metavar='p1', type=float, default=-0.00000102, required=False,
                              help=argparse.SUPPRESS)
    pbseq_parser.add_argument("--sp2", metavar='p2', type=float, default=0.00004281, required=False,
                              help=argparse.SUPPRESS)

    args_config, _ = parser.parse_known_args()

    if args_config.config:
        config_values = load_config_file(args_config.config)
        parser.set_defaults(**config_values)

        # Set arguments provided via config as no longer required,
        # this way the parser doesn't throw an error when parsing command line args.
        for action in parser._actions:
            if action.dest in config_values:
                action.required = False

    args = parser.parse_args()
    validate_args(args)
    args = infer_defaults(args)

    return args


def validate_args(args):

    # This branch is here to allow for input-seq to be overridden as "None" if it's provided in a configuration
    # file but the user desires to provide a input-pop file.
    if args.input_seq in ("None", "none", ""):
        args.input_seq = None

    if args.reads is None:
        raise ValueError(f'--reads is a required input.')

    if args.rtp is None:
        raise ValueError(f'--rtp is a required input.')

    if args.ncdna is None:
        raise ValueError(f'--ncdna is a required input.')

    if args.input_seq is None and args.input_pop is None and args.random is None:
        args_debug_str = print_selected_args(args, ['input_seq', 'input_pop', 'random'])
        raise ValueError(f'None of --input-seq, --input-pop, or --random provided. Please check inputs.\n'
                         f'{args_debug_str}')
    if args.input_pop and args.input_seq:
        args_debug_str = print_selected_args(args, ['input_seq', 'input_pop'])
        raise ValueError('Conflicting inputs --input-seq and --input-pop cannot be used together.\n'
                         f'{args_debug_str}')
    if (args.input_pop or args.input_seq) and args.random:
        args_debug_str = print_selected_args(args, ['input_seq', 'input_pop', 'random'])
        raise ValueError('--random cannot be used with a provided user sequence.\n'
                         f'{args_debug_str}')

    if args.dirichlet:
        if args.input_seq is None and args.random is None:
            args_debug_str = print_selected_args(args, ['dirichlet', 'input_seq', 'random'])
            raise ValueError('--dirichlet was provided without an input sequence configuration.\n'
                             f'{args_debug_str}')
        if args.input_pop:
            args_debug_str = print_selected_args(args, ['dirichlet', 'input_pop'])
            raise ValueError('--dirichlet is not compatible with --input-pop; it requires --input-seq or --random.\n'
                             f'{args_debug_str}')

        if args.clonal:
            args_debug_str = print_selected_args(args, ['dirichlet', 'clonal'])
            raise ValueError('--dirichlet is not compatible with --clonal.\n'
                             f'{args_debug_str}')

        if args.nhaps is None:
            args_debug_str = print_selected_args(args, ['dirichlet', 'nhaps'])
            raise ValueError('--dirichlet requires --nhaps.\n'
                             f'{args_debug_str}')

    if args.clonal:
        if args.input_pop is not None:
            args_debug_str = print_selected_args(args, ['input_pop', 'clonal'])
            raise ValueError('--input-pop is not compatible with --clonal.\n'
                             f'{args_debug_str}')

        if args.ncdna is None:
            args_debug_str = print_selected_args(args, ['ncdna', 'clonal'])
            raise ValueError('--clonal requires --ncdna to be provided.\n'
                             f'{args_debug_str}')

        if args.nhaps is not None:
            args_debug_str = print_selected_args(args, ['nhaps', 'clonal'])
            raise ValueError('--clonal cannot be used with --nhaps to be provided.\n'
                             f'{args_debug_str}')

    if args.error_rate_rt >= 1.0 or args.error_rate_rt < 0.0:
        args_debug_str = print_selected_args(args, ['error_rate_rt'])
        raise ValueError('--error_rate_rt must be >=0 and <1.\n'
                         f'{args_debug_str}')

    if args.cycles < 1:
        args_debug_str = print_selected_args(args, ['cycles'])
        raise ValueError('--cycles must be >=1.\n'
                         f'{args_debug_str}')

    if args.input_pop is None:
        if args.ncdna is None:
            args_debug_str = print_selected_args(args, ['input_pop', 'ncdna'])
            raise ValueError('--ncdna must be provided unless --input-pop is given.\n'
                             f'{args_debug_str}')

    if args.input_pop is not None:
        if args.dirichlet is not None:
            args_debug_str = print_selected_args(args, ['input_pop', 'dirichlet'])
            raise ValueError('--input_pop cannot be provided with --dirichlet.\n'
                             f'{args_debug_str}')

        if args.nhaps is not None:
            args_debug_str = print_selected_args(args, ['input_pop', 'nhaps'])
            raise ValueError('--input_pop cannot be provided with --nhaps.\n'
                             f'{args_debug_str}')

    if args.fp is not None and args.sss is True:
        args_debug_str = print_selected_args(args, ['fp', 'sss'])
        raise ValueError('--fp cannot be provided with --sss. Provide forward primers via --sss-fp.\n'
                         f'{args_debug_str}')

    if args.error_rate_sss is None:
        args.error_rate_sss = args.error_rate_pcr

    if args.recomb_rate_sss is None:
        args.recomb_rate_sss = args.recomb_rate_pcr

    if args.movie_time < 10 or args.movie_time > 40:
        args_debug_str = print_selected_args(args, ['movie_time'])
        raise ValueError('--movie-time must be between 10 and 40 hours.\n'
                         f'{args_debug_str}')


def infer_defaults(args):
    if args.name is None:
        if args.input_seq:
            args.name = Path(args.input_seq).stem
        if args.input_pop:
            args.name = Path(args.input_pop).stem
        if args.name is None:
            args.name = 'umisim'
        args.name = datetime.datetime.now().strftime(f'{args.name}-%y%m%d-%H%M%S')

    hex_id = "{:03x}".format(random.randint(0, 0xFFF)) + "{:01x}".format(random.randint(0xA, 0xF))
    args.hex_id = hex_id

    if args.output is None:
        args.output = Path(args.name).absolute()
    else:
        args.output = Path(args.output).absolute()

    return args


def load_config_file(fp_config):
    if fp_config.endswith(('.yaml', '.yml')):
        with open(fp_config, 'r') as f:
            return yaml.safe_load(f)
    else:
        raise ValueError(f'Configuration file must have .yaml/.yml extension.')


def print_selected_args(args, selected_args):
    msg = ''
    for selected_arg in selected_args:
        parg = selected_arg.replace('_', '-')
        msg += f'\n\t--{parg}={args.__dict__[selected_arg]}'
    return msg

# if __name__ == '__main__':
#     multiprocessing.freeze_support()
