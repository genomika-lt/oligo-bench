"""Describes function for adding id to csv files"""


import logging

import pandas as pd

from snakemake.script import snakemake


# logging
logger = logging.getLogger(__name__)
logging.basicConfig(filename=snakemake.log[0],
                    filemode='w',
                    encoding='utf-8',
                    level=logging.INFO)


def add_id_to_csv(input_file, column_id, samples, output_file):
    """
    Adds id to csv file.
    :param str input_file: Input file to modify
    :param column_id: id of column
    :param samples: value to insert
    :param str output_file: Output file to save data
    :rtype: None
    """

    df = pd.read_csv(input_file)
    df.insert(0, column_id, samples)
    df.to_csv(output_file, index=False)


try:
    add_id_to_csv(input_file=snakemake.input[0], column_id=snakemake.params.id,
              samples=snakemake.params.samples, output_file=snakemake.output[0])
except Exception as e:
    logger.exception(e)
    raise e
