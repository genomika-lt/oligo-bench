"""Converts TXT to FASTA reference file"""

from snakemake.script import snakemake

# pylint: disable=import-error
from scripts.utils import snakemake_file_logger


@snakemake_file_logger
def convert_txt_to_fasta(input_txt, output_fasta):
    """
    Convert a text file to FASTA format with zero-padded sequence numbers.
    :param input_txt: Path to the input .txt file
    :param output_fasta: Path to the output .fa file
    """
    with open(input_txt, 'r') as txt_file:
        with  open(output_fasta, 'w') as fasta_file:
            for line_number, line in enumerate(txt_file, start=1):
                sequence = line.strip()
                if sequence:
                    fasta_file.write(f">{line_number:06d}\n")
                    fasta_file.write(f"{sequence}\n")


convert_txt_to_fasta(input_txt=snakemake.params['input_txt'], output_fasta=snakemake.output[0])