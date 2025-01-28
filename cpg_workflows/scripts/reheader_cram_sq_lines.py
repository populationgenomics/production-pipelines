import os
import subprocess as sb
import sys
from argparse import ArgumentParser


def cli_main():
    parser = ArgumentParser(description='CLI for the CRAM @SQ line reheadering script')
    parser.add_argument('--cram_in', help='Path to a localised CRAM, this will be modified in place', required=True)
    parser.add_argument('--sq_file', help='Path to a localised file containing the new @SQ header lines', required=True)
    args = parser.parse_args()

    reheader_cram_sq_lines(
        file_in=args.cram_in,
        sq_file=args.sq_file,
    )

def reheader_cram_sq_lines(
    input_cram: str, sq_file: str,
):
    """
    Reads the CRAM header lines and replaces each @SQ line with the corresponding lines from the SQ file.
    Leaves all other header lines unchanged.

    args:
        input_cram: str, path to the input CRAM file
        output_cram: str, path to the output CRAM file
        sq_file: str, path to the SQ header lines file
    """
    header_lines, sq_index = read_headers(input_cram)
    if header_lines is None:
        print(f'Error: Could not read the headers from {input_cram}', file=sys.stderr)
        sys.exit(1)

    sq_lines = read_sq_lines(sq_file)
    if sq_lines is None:
        print(f'Error: Could not read the @SQ lines from {sq_file}', file=sys.stderr)
        sys.exit(1)

    new_header_lines = add_sq_lines_to_header(header_lines, sq_index, sq_lines)
    new_header_file = os.path.join(os.environ['BATCH_TMPDIR'], 'tmp.hdr')
    with open(new_header_file, 'w') as new_header:
        for line in new_header_lines:
            print(line, file=new_header)


    run = sb.run(
        ['samtools', 'reheader', '--no-PG', '--in-place', new_header_file, input_cram],
        check=True,
    )
    if run.returncode != 0:
        print(f'Error: Could not reheader the CRAM file: {run.stderr}', file=sys.stderr)
        sys.exit(1)
    else:
        print(f'Successfully reheadered the CRAM file: {input_cram}')

    run_index = sb.run(
        ['samtools', 'index', input_cram],
        check=True,
    )
    if run_index.returncode != 0:
        print(f'Error: Could not index the CRAM file: {run_index.stderr}', file=sys.stderr)
        sys.exit(1)
    else:
        print(f'Successfully indexed the CRAM file: {input_cram}.crai')


def add_sq_lines_to_header(cram_header: list[str], sq_index: int, sq_lines: list[str]):
    """
    Replaces the @SQ lines in the CRAM header with the @SQ lines from the SQ file.

    args:
        cram_header: list[str], the CRAM header lines
        sq_index: int, the index at which the new @SQ lines should be inserted
        sq_lines: list[str], the new @SQ lines
    """
    return cram_header[:sq_index] + sq_lines + cram_header[sq_index + len(sq_lines):]


def read_headers(path: str):
    """Reads the file's headers and returns all lines that are not @SQ lines, and the index at which the new @SQ lines should be inserted."""
    run = sb.run(['samtools', 'head', path], capture_output=True, text=True)
    if run.returncode != 0:
        return None  # Probably file not found

    sq_index = None
    all_header_lines = []
    for line in run.stdout.splitlines():
        if not line.startswith('@SQ'):
            all_header_lines.append(line)
        else:
            sq_index = sq_index or len(all_header_lines)

    return all_header_lines, sq_index


def read_sq_lines(path: str):
    """Reads the @SQ lines from the SQ file."""
    with open(path) as f:
        return [line for line in f if line.startswith('@SQ')]