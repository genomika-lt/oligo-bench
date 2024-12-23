"""Creates mapping table summary"""


from pandas import DataFrame
from snakemake.script import snakemake

# pylint: disable=import-error
from scripts.utils import (parse_sam_bam_file,
                           snakemake_file_logger,
                           integer_to_human_readable)


@snakemake_file_logger
def payload_errors_number(sam_bam_files, output_file):
    """
    Counts errors in aligned reads from bam files and gives ratio.
    :param list[str] sam_bam_files: List with paths to SAM/BAM files.
    :param str output_file: Output file to save data.
    :rtype: None
    """


    data = {'Sample': [],
            'Payload Errors': []}

    for file in sam_bam_files[:len(sam_bam_files) // 2]:
        total_error = 0
        total_reads = 0
        for read in parse_sam_bam_file(file):
            total_reads += 1
            total_error += read.get_tag('NM') / len(read.query_qualities)

        sample = file.split('/')[-1].split('.')[0][9:]
        data['Sample'].append(sample)
        data['Payload Errors'].append(total_error / total_reads)

    data_frame = DataFrame(data)
    data_frame.to_csv(output_file, index=False)


payload_errors_number(sam_bam_files=snakemake.input, output_file=snakemake.output[0])
