"""File collecting most useful functions as utils."""


import logging
from typing import Callable
from math import floor, log10

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


def round_to_x_significant(number, x):
    """

    :param number:
    :param x:
    :return:
    """
    return round(number, x - 1 - floor(log10(abs(number))))


def from_phred_to_char(value: int|list[int]) -> str:
    if isinstance(value, int):
        return chr(value + 33)

    return ''.join([chr(i + 33) for i in value])


def from_char_to_phred(value: str|list[str]) -> int|list[int]:
    if len(value) == 1:
        return ord(value[0]) - 33

    return [ord(i) - 33 for i in value]
