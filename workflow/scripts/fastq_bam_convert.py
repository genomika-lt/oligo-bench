"""Converts FASTQ to BAM file"""

import numpy as np
import pysam
import gzip
import shutil

from snakemake.script import snakemake

# pylint: disable=import-error
from scripts.utils import snakemake_file_logger


def calculate_qscore(qstring: str):
    """
    Calculates a mean_qscore from a qstring
    :param qstring used for getting mean qscore value
    :return: mean_qscore
    """
    qs = (np.array(qstring, 'c').view(np.uint8) - 33)
    mean_err = np.exp(qs * (-np.log(10) / 10.)).mean()
    mean_qscore = -10 * np.log10(max(mean_err, 1e-4))
    return mean_qscore


def convert_fastq_to_bam(fastq_files: list[str], bam_file: str):
    """
    Converts FASTQ file to BAM file
    :param fastq_files: FASTQ file used for input conversion
    :param bam_file: BAM file used for result of conversion
    :return: None
    """
    with pysam.AlignmentFile(bam_file, "wb", header={"HD": {"VN": "1.6"}}) as bam:
        for file in fastq_files:
            with open(file, "r") as fq:
                while True:
                    seq_id = fq.readline().strip()
                    if not seq_id:
                        break

                    sequence_meta = seq_id.split(" ")
                    metadata = {"seq_id": sequence_meta[0]}
                    for part in sequence_meta[1:]:
                        key, value = part.split("=")
                        metadata[key] = value

                    sequence = fq.readline().strip()
                    fq.readline().strip()
                    quality = fq.readline().strip()

                    read = pysam.AlignedSegment()
                    read.query_name = metadata["seq_id"]
                    read.query_sequence = sequence
                    read.flag = 4
                    read.reference_id = -1
                    read.reference_start = -1
                    read.mapping_quality = 0
                    read.cigar = tuple()
                    read.query_qualities = pysam.qualitystring_to_array(quality)

                    read.set_tag("RG", metadata.get("protocol_group_id"))
                    read.set_tag("QS", calculate_qscore(quality))
                    read.set_tag("ST", metadata.get("start_time"))
                    read.set_tag("SM", metadata.get("sample_id"))
                    read.set_tag("BC", metadata.get("basecall_model_version_id"))
                    read.set_tag("PU", metadata.get("flow_cell_id"))
                    read.set_tag("CO", metadata.get("parent_read_id"))
                    read.set_tag("CR", metadata.get("runid"))
                    read.set_tag("CN", metadata.get("ch"))

                    bam.write(read)


@snakemake_file_logger
def provide_bam(fastq_files: list[str],
                    bam_files: list[str],
                    bam_file: str):
    if bam_files:
        return merge_bam_files(bam_files, bam_file)
    if fastq_files:
        fastq_files = unzip_fastq_files(fastq_files)
        return convert_fastq_to_bam(fastq_files, bam_file)
    sample_id = bam_file.split("/")[-1].split(".")[0]
    raise FileNotFoundError(f"No FASTQ or BAM files found in {sample_id}")


def unzip_fastq_files(fastq_files: list[str]) -> list[str]:
    """
    Ensures all FASTQ files are uncompressed. If a file is .gz, unzip it.
    :param fastq_files: list of FASTQ files
    :return: list of uncompressed FASTQ files
    """
    unzipped_files = []
    for file in fastq_files:
        if file.endswith(".gz"):
            unzipped_file = file[:-3]
            with gzip.open(file, 'rb') as f_in:
                with open(unzipped_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            unzipped_files.append(unzipped_file)
        else:
            unzipped_files.append(file)
    return unzipped_files


def merge_bam_files(bam_files: list[str], bam_file: str):
    """
    Merges list of BAM files to one BAM file
    :param list[str] bam_files: list of BAM files
    :param str bam_file: Output file to save data
    :rtype: None
    """
    sam_file = pysam.AlignmentFile(bam_files[0], "rb",check_sq=False)

    with pysam.AlignmentFile(bam_file, "wb",template=sam_file) as out_bam:
        for bf in bam_files:
            with pysam.AlignmentFile(bf, "rb",check_sq=False) as in_bam:
                for read in in_bam:
                    out_bam.write(read)


provide_bam(fastq_files=snakemake.params['fastq_files'],
                bam_files=snakemake.params['bam_files'],
                bam_file=snakemake.output[0])