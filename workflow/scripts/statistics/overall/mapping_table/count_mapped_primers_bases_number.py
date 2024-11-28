"""Creates mapping table summary"""


from pandas import DataFrame
from snakemake.script import snakemake

# pylint: disable=import-error
from scripts.utils import (parse_sam_bam_file,
                           snakemake_file_logger)


@snakemake_file_logger
def count_mapped_primers_reads_number(sam_bam_files, output_file):
    """
    Counts total number of mapped bases from sam files and saves as csv.
    :param list[str] sam_bam_files: List with paths to SAM/BAM files.
    :param str output_file: Output file to save data.
    :rtype: None
    """


    data = {'Sample': [],
            'Mapped Bases': []}

    for file in sam_bam_files:
        mapped_bases = 0

        for read in parse_sam_bam_file(file):
            mapped_bases += len(read.query_qualities)

        sample = file.split('/')[-1].split('.')[0][16:]
        data['Sample'].append(sample)
        data['Mapped Bases'].append(mapped_bases)

    data_frame = DataFrame(data)
    data_frame.to_csv(output_file, index=False)


count_mapped_primers_reads_number(sam_bam_files=snakemake.input, output_file=snakemake.output[0])
