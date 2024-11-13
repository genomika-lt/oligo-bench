""" Counts number of passed bases inside SAM/BAM files """


from pandas import DataFrame
from snakemake.script import snakemake

# pylint: disable=import-error
from scripts.utils import (parse_sam_bam_file,
                           snakemake_file_logger)

# pylint: disable=too-many-locals
@snakemake_file_logger
def count_passed_bases_number(sam_bam_files, output_file):
    """
    Counts number of passed bases inside SAM/BAM files and saves as csv file.
    :param list[str] sam_bam_files: List with paths to SAM/BAM files.
    :param str output_file: Output file to save data.
    :rtype: None
    """


    data = {'Sample': [],
            'Passed Bases': []}

    for file in sam_bam_files:
        bases_counter = 0
        for read in parse_sam_bam_file(file):
            bases_counter += len(read.query_sequence)

        sample = file.split('/')[-1].split('.')[0]
        data['Sample'].append(sample)
        data['Passed Bases'].append(bases_counter)

    data_frame = DataFrame(data)
    data_frame.to_csv(output_file, index=False)


count_passed_bases_number(sam_bam_files=snakemake.input, output_file=snakemake.output[0])
