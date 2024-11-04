"""Saves primers in fasta file for last aligner"""


from snakemake.script import snakemake

# pylint: disable=import-error
from scripts.utils import snakemake_file_logger


@snakemake_file_logger
def extract_primers(forward_primer: str, reverse_primer: str, output_file):
    """
    Plots gc distribution over time in samples and saves to html
    :param str forward_primer: Forward primer sequence as string.
    :param str reverse_primer: Reverse primer sequence as string.
    :param str output_file: Output file to save primers to as fasta file.
    :rtype: None
    """
    complemented_dictionary = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_primer_complemented = ''.join([complemented_dictionary[i]
                                           for i in reverse_primer][::-1])

    with open(output_file, "w", encoding='utf-8') as g:
        g.write(">Forward5\n")
        g.write(forward_primer)

        g.write("\n")

        g.write(">Reverse3\n")
        g.write(reverse_primer_complemented)


extract_primers(forward_primer=snakemake.params['forward_primer'],
                reverse_primer=snakemake.params['reverse_primer'],
                output_file=snakemake.output[0])
