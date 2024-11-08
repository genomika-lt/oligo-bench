"""File collecting most useful functions as utils."""


import logging
from collections.abc import Iterable
from math import floor, log10

import pysam
from snakemake.script import snakemake


def snakemake_file_logger(func):
    """
    Decorator for logging snakemake functions into files.
    :param Callable func: Snakemake script function to wrap.
    :return: Wrapped function.
    :rtype: Callable
    """


    logging.basicConfig(filename=snakemake.log[0],
                        filemode='w',
                        encoding='utf-8',
                        level=logging.INFO)


    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            logging.exception(e)
            raise e


    return wrapper


def parse_sam_records(records):
    """
    Parses bam records.
    :param list[list[str]] records: Records from BAM file.
    :return: Parsed records.
    """


    parsed_records = []
    for record in records:
        meta = {}
        for meta_data in record[11:]:
            key_index = meta_data.index(':')
            type_index = meta_data.index(':', key_index + 1)
            meta_key, meta_type, meta_value = (meta_data[:key_index],
                                               meta_data[key_index + 1:type_index],
                                               meta_data[type_index + 1:])
            if meta_type == 'i':
                meta[meta_key] = int(meta_value)
            elif meta_type == 'f':
                meta[meta_key] = float(meta_value)
            elif meta_type == 'Z':
                meta[meta_key] = meta_value
            else:
                meta[meta_key] = meta_value

        parsed_record = record[:11] + [meta]
        parsed_records.append(parsed_record)

    return parsed_records

# pylint: disable=no-member
def parse_bam_file(path_to_file: str) -> pysam.AlignedSegment:
    """
    Generator that parses bam file and yields read by read.
    :param path_to_file: Path to BAM file.
    :return: Yields read from BAM file.
    """


    reads: pysam.AlignmentFile = pysam.AlignmentFile(path_to_file, 'rb', check_sq=False)
    yield from reads


def round_to_x_significant(number: int|float, x: int):
    """
    Rounds number to x significant digits.
    :param number: Number to round.
    :param x: number of digits.
    :return: Rounded number.
    """


    return round(number, x - 1 - floor(log10(abs(number))))


def from_phred_to_char(value: int|Iterable[int]) -> str:
    """
    Converts integer value to character or list of integers to list of characters.
    :param value: Phred value or list of phred values.
    :return: Character of list of characters.
    """


    if isinstance(value, int):
        return chr(value + 33)

    return ''.join([chr(i + 33) for i in value])


def from_char_to_phred(value: str|Iterable[str]) -> int|list[int]:
    """
    Converts character to phred value or list of characters to list of values.
    :param value: Character of list of characters.
    :return: Phred value or list of phred values.
    """


    if len(value) == 1:
        return ord(value[0]) - 33

    return [ord(i) - 33 for i in value]


def integer_to_human_readable(number: int,
                             alphabet: list[str] | tuple[str] = ('', 'K', 'M', 'B', 'T')) -> str:
    """
    Converts integer value to human-readable string with using user specified alphabet.
    :param number: Number to convert.
    :param alphabet: Alphabet to use.
    :return: Converted, readable value.
    """


    string_number = str(number)
    number_of_digits = len(string_number)
    suffix_index = (number_of_digits - 1) // 3

    if suffix_index >= len(alphabet):
        raise IndexError(f'Too small alphabet: {len(alphabet)},'
                         f' was referring to element number {suffix_index}')

    return string_number[:len(string_number) - suffix_index * 3] + alphabet[suffix_index]
