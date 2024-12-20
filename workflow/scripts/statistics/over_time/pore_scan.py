"""Plots pore activity based on minknow output csv file"""


import os
import math
import json

import plotly.graph_objects as go
from plotly.subplots import make_subplots

from snakemake.script import snakemake

# pylint: disable=import-error
from scripts.utils import snakemake_file_logger


# pylint: disable=too-many-locals
@snakemake_file_logger
def pore_scan(path_to_samples, output_file):
    """
    Plots pore activity based on minknow output csv file
    :param list[str] path_to_samples: List of paths to samples
    :param str output_file: Output file to save data
    :rtype: None
    """

    json_reports = []
    for path in path_to_samples:
        for file in os.listdir(path):
            if file.endswith(".json"):
                json_reports.append(os.path.join(path, file))

    samples_scans = []
    number_of_pores = []
    for report in json_reports:
        with open(report, 'r', encoding='utf-8') as f:
            json_data = json.loads(f.read())

        sample_pore_scans: list[dict] = json_data['acquisitions'][-1]['acquisition_run_info'][
            'bream_info']['mux_scan_results']
        flow_cell_number_of_pores = (json_data['protocol_run_info']['flow_cell']['channel_count'] *
                                     json_data['protocol_run_info']['flow_cell']['wells_per_channel'])

        number_of_pores.append(flow_cell_number_of_pores)
        samples_scans.append(sample_pore_scans)

    max_number_of_scans = 0
    for sample in samples_scans:
        max_number_of_scans = max(max_number_of_scans, len(sample))

    pore_types = set()
    for sample in samples_scans:
        pore_types.update(set(sample[0]['counts'].keys()))

    sample_names = [path.split('/')[-2] for path in path_to_samples]
    mapped_number_of_pores = {}
    for i in range(len(sample_names)):
        mapped_number_of_pores[sample_names[i]] = number_of_pores[i]

    pore_scans = [{pore: {sample_names[sample]: (samples_scans[sample][scan]['counts'][pore]
                                            if len(samples_scans[sample]) > scan else 0)
                          for sample in range(len(sample_names))}
                   for pore in pore_types}
                  for scan in range(max_number_of_scans)]

    possible_number_of_columns = [i for i in range(1, 48 + 1) if 49 % i == 0][::-1]
    max_number_columns = 48 / len(possible_number_of_columns)
    for possible_value in possible_number_of_columns:
        if max_number_columns >= possible_value:
            max_number_columns = possible_value
            break

    actual_number_of_columns = len(pore_scans)
    if len(pore_scans) >= max_number_columns:
        actual_number_of_columns = max_number_columns

    actual_number_of_rows = math.ceil(len(pore_scans) / max_number_columns)

    pore_types_metadata = json_data['acquisitions'][-1]['acquisition_run_info']['bream_info'][
        'mux_scan_metadata']['category_groups']
    pore_type_color = {pore_type['name']: pore_type['style']['colour']
                       for pore_type in pore_types_metadata}
    pore_type_color['multiple'] = 'ff000f'
    pore_type_color['other'] = '000000'

    pore_types_ordered = ['single_pore',
                          'reserved_pore',
                          'unavailable',
                          'saturated',
                          'zero',
                          'multiple',
                          'other']

    if len(pore_types.union(pore_types_ordered)) != len(pore_types_ordered):
        raise ValueError('Pore types between current version of QC and minknow do not match')
    figure = make_subplots(rows=actual_number_of_rows,
                           cols=actual_number_of_columns,
                           shared_yaxes=True)

    for scan_index, scan_value in enumerate(pore_scans):
        for pore_type in pore_types_ordered:
            figure.add_trace(go.Bar(x=list(scan_value[pore_type].keys()),
                                    y=[scan_value[pore_type][sample] / mapped_number_of_pores[sample]
                                       for sample in scan_value[pore_type]],
                                    name=pore_type,
                                    legendgroup=pore_type,
                                    hovertext=f'{pore_type}',
                                    showlegend=scan_index == 0,
                                    marker_color='#' + pore_type_color[pore_type].lower()),
                             row=scan_index // actual_number_of_columns + 1,
                             col=scan_index % actual_number_of_columns + 1)

    figure.update_layout(barmode='stack')

    figure.update_layout(title="Pore Scan results",
                         legend_title="Pore types")

    with open(output_file, 'w', encoding='utf-8') as g:
        g.write(figure.to_html(full_html=False, include_plotlyjs='cdn'))


pore_scan(path_to_samples=snakemake.input, output_file=snakemake.output[0])
