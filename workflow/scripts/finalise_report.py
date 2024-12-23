"""Counts total number of reads"""


from io import TextIOWrapper

from snakemake.script import snakemake

# pylint: disable=import-error
from scripts.utils import snakemake_file_logger


def write_summary_to_file(file: TextIOWrapper) -> None:
    params = snakemake.config['basecalling']
    file.write("<details>")
    file.write('<summary>Information</summary>')
    file.write("<p>QC Version: v1.2</p>")
    file.write("<p>Dorado: 0.9.0</p>")
    file.write(f"<p>Dorado model: {params['dorado_model']['value']}</p>")
    file.write(f"<p>Minimum basecalling quality: {params['minimum_quality']['value']}</p>")
    file.write(f"<p>Allowed primer errors: {params['allowed_primer_errors']['value']}</p>")
    file.write("</details>\n")


@snakemake_file_logger
def finalise_report(reports, output_file):
    """
    Counts total number of reads across samples
    :param str reports: Folder with fastq reads
    :param str output_file: Output file to save data
    :rtype: None
    """

    with open(output_file, 'w', encoding='utf-8') as g:
        write_summary_to_file(g)
        for report in reports:
            with open(report, 'r', encoding='utf-8') as f:
                g.write(f.read())


finalise_report(reports=snakemake.input, output_file=snakemake.output[0])
