""" Counts cumulative number of bases inside BAM files """

import os
import sys
import json
from datetime import datetime, timedelta
import plotly.graph_objects as go
import pysam
from snakemake.script import snakemake
from scripts.utils import snakemake_file_logger

# pylint: disable=too-many-locals
@snakemake_file_logger
def cumulative_bases_plot(path_to_samples, sorted_bam_file, output_file):
    """
    Generates a cumulative number of bases plot based
    on sample data and flow cell specifications.

    This function reads JSON report files to gather data
    on theoretical and practical estimates based on flow cell type
    and acquisition data. It uses information
    from the BAM file to calculate
    the cumulative number of bases over time.
    The plot compares the actual cumulative bases (real
    values) with the theoretical and practical estimates.

    :param str path_to_samples: A path to samples to be analyzed.
    :param str sorted_bam_file: A path to sorted BAM file.
    :param str output_file: The path to the output file where
     the generated plot will be saved as an HTML file.

    Raises:
        ValueError: If an unknown flow cell product code is encountered in the JSON data.
    """

    flow_cell_data = {
        "flg": {"gb": 2.8,"channels": 126},
        "minion": {"gb": 50,"channels": 512},
        "promethion": {"gb": 290,"channels": 3000},
    }

    json_reports = []
    for path in path_to_samples:
        for file in os.listdir(path):
            if file.endswith(".json"):
                json_reports.append(os.path.join(path, file))

    for report in json_reports:
        with open(report, 'r', encoding='utf-8') as f:
            json_data = json.loads(f.read())

    fig = go.Figure()
    data = count_bases(str(sorted_bam_file))

    time_points = list(data.keys())
    product_code = (
        json_data['protocol_run_info']['flow_cell']
        ['user_specified_product_code'].lower())

    match product_code:
        case code if "flg" in code:
            flow_cell_type = "flg"
        case code if "min" in code:
            flow_cell_type = "minion"
        case code if "prom" in code:
            flow_cell_type = "promethion"
        case _:
            raise ValueError("Unknown product code: Unable to calculate coefficient")

    coefficient = flow_cell_data[flow_cell_type]["gb"] / flow_cell_data[flow_cell_type]["channels"]
    flow_cell_gb = flow_cell_data[flow_cell_type]["gb"]

    theoretical_values =[flow_cell_gb*1e9] * len(time_points)
    practical_values = [(json_data['acquisitions'][-1]['acquisition_run_info']['bream_info']
                         ['mux_scan_results'][0]['counts']['single_pore']+
                        json_data['acquisitions'][-1]['acquisition_run_info']['bream_info']
                        ['mux_scan_results'][0]['counts']
                        ['reserved_pore'])*coefficient*1e9] * len(time_points)
    real_values = []
    real_values_total = 0
    for value in data.values():
        real_values_total += value
        real_values.append(real_values_total)

    fig.add_trace(go.Scatter(
        x=time_points,
        y=theoretical_values,
        mode='lines',
        name=f'{flow_cell_type.upper()} ({flow_cell_data[flow_cell_type]["gb"]} Gb) t.n.b',
        line={"color": 'red', "dash": 'dot'}
    ))

    fig.add_trace(go.Scatter(
        x=time_points,
        y= practical_values,
        mode='lines',
        name=f'{flow_cell_type.upper()} p.n.b',
        line={"color": 'blue', "dash": 'dash'}
    ))

    fig.add_trace(go.Scatter(
        x=time_points,
        y=real_values,
        mode='lines',
        name=f'{flow_cell_type.upper()} r.n.b',
        line={"color": 'green', "dash": 'solid'}
    ))

    fig.update_layout(
        title="Cumulative number of bases",
        xaxis_title="Time(sec)",
        yaxis_title="Number of Bases (Gb)",
        xaxis={"tickformat": '%H:%M', "title_standoff": 25},
        yaxis={"title_standoff": 25, "type": 'log', "tickformat": '.2e'},
        margin={"l": 50, "r": 50, "t": 50, "b": 50}
    )

    with open(output_file, 'w', encoding='utf-8') as g:
        g.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))


def count_bases(bam_file, batch_size: int = 1073741824):
    """
    Reads a BAM file, processes it in batches, and counts the total number of bases in each batch.

    This function iterates over the sorted BAM file to group reads
    into batches. The batch is formed by reads that
    have the same timestamp and whose total size does not exceed
    the specified `batch_size`. Once a batch is formed,
    it calculates the total number of bases in the batch and
    stores the result with the timestamp of the first read in that batch.

    :param str bam_file: Path to the input BAM file to be processed.
    :param int batch_size: Maximum size of each batch (in bytes). Default is 1073741824 bytes.

    Returns:
        dict: A dictionary where keys are timestamps and values are the total number of bases
              for each corresponding batch.

    Raises:
        ValueError: If there is an issue with reading the BAM file or extracting the required tags.
    """
    num_bases={}
    with pysam.AlignmentFile(bam_file, "rb", check_sq=False) as bam:
        read = next(bam)
        batch = [read]
        current_batch_size = sys.getsizeof(read)
        try:
            while True:
                read = next(bam)
                read_size = sys.getsizeof(read)
                first_st = datetime.fromisoformat(batch[0].get_tag('st'))
                current_st = datetime.fromisoformat(read.get_tag('st'))
                if (abs(first_st-current_st) <= timedelta(seconds=60) and
                        current_batch_size + read_size <= batch_size):
                    batch.append(read)
                    current_batch_size += read_size
                else:
                    total_sum = sum(len(read.query_sequence) for read in batch)
                    time_stamp = datetime.fromisoformat(batch[0].get_tag('st'))
                    num_bases[time_stamp] = total_sum
                    batch = [read]
                    current_batch_size = read_size
        except StopIteration:
            pass

    return num_bases

cumulative_bases_plot(path_to_samples=snakemake.input["path_to_samples"],
                      sorted_bam_file=snakemake.input["bam_file"],
                      output_file=snakemake.output[0])
