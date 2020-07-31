import os
import subprocess


def cdhit(
        input_fasta: str,
        percent_identity: float,
        length_difference_cutoff: float = None,
        min_alignment_coverage: float = None,
        percent_identity_suffix: bool = False,
        output_dir: str = None
) -> None:
    """Runs the cd-hit program.

    Args:
        input_fasta:
        output_dir: (Optional) If provided, writes cd-hit files to this directory

    Returns:
        A bool indicate True if program was run with non zero error code.

    Raises:
        ValueError when...
        CalledProcessError when cd-hit completes with an error.
    """

    if not os.path.exists(input_fasta):
        raise FileNotFoundError(f'Fasta file {input_fasta} was not found.')

    if output_dir is not None:
        if not os.path.isdir(output_dir):
            raise NotADirectoryError(f'{output_dir} is not a directory.')
        output_prefix = os.path.join(output_dir, os.path.splitext(os.path.basename(input_fasta))[0])
    else:
        output_prefix = os.path.splitext(input_fasta)[0]

    if percent_identity_suffix:
        output_prefix = output_prefix + str(percent_identity)
    # from pdb import set_trace; set_trace()
    params = [
        'cd-hit',
        '-i', input_fasta,
        '-o', output_prefix,
        '-c', str(percent_identity),
        '-n', str(_cdhit_word_size(percent_identity)),
        '-d', '0',
        '-sc', '1',
    ]

    if length_difference_cutoff is not None:
        params.extend(['-s', str(length_difference_cutoff)])

    if min_alignment_coverage is not None:
        params.extend([
            '-aL', str(min_alignment_coverage),
            '-aS', str(min_alignment_coverage),
        ])

    subprocess.check_call(params)


def _cdhit_word_size(percent_id):
    if percent_id < .4:
        raise ValueError("No word size or percent identity < 0.4")
    elif 0.4 <= percent_id < .5:
        return 2
    elif 0.5 <= percent_id < .6:
        return 3
    elif 0.6 <= percent_id < .7:
        return 4
    elif 0.7 <= percent_id <= 1.0:
        return 5
    else:
        raise ValueError("Percent identity for cd-hit must 1.0 or less.")
