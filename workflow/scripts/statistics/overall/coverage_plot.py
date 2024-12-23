"""Plots coverage for each oligo"""

import plotly.graph_objects as go

from snakemake.script import snakemake

# pylint: disable=import-error
from scripts.utils import snakemake_file_logger, parse_sam_bam_file


@snakemake_file_logger
def coverage_plot(bam_ref_files, output_file):
    """
    Plots sorted coverage for each oligo.
    :param list[str] bam_ref_files: list with bam and reference files.
    :param str output_file: Output file to save data.
    :rtype: None
    """

    # pylint: disable=duplicate-code
    figure = go.Figure()

    number_of_bam_files = len(bam_ref_files) // 2
    for bam_file, ref_file in zip(bam_ref_files[:number_of_bam_files],
                                  bam_ref_files[number_of_bam_files:]):
        coverage = {}
        with open(ref_file, encoding='utf-8') as ref:
            for line in ref.readlines():
                if line.startswith('>'):
                    coverage[line.strip()[1:]] = 0

        for read in parse_sam_bam_file(bam_file):
            coverage[read.reference_name] += 1

        data = []

        for key, value in coverage.items():
            data.append((key, value))

        data.sort(key=lambda x: x[1], reverse=True)

        y = [i[1] for i in data]

        stretched_y: list[float] = [0.0 for _ in range(100)]

        for i in range(100):
            left_index = (i * len(data)) // 100
            right_index = ((i + 1) * len(data)) // 100
            if left_index == right_index:
                stretched_y[i] = 0

            stretched_y[i] = sum(y[left_index:right_index:right_index]) / (right_index - left_index)

        max_coverage = max(stretched_y)
        if max_coverage == 0:
            max_coverage = 1

        normalized_y = [i / max_coverage for i in stretched_y]

        figure.add_trace(go.Scatter(y=normalized_y,
                                    x=[i for i in range(100)],
                                    name=bam_file.split('/')[-1][7:-4]))
        data.sort(key=lambda x: x[1], reverse=True)

    figure.update_layout(title="Passed Aligned Reads Coverage",
                         xaxis_title="Oligo",
                         yaxis_title="Coverage",
                         legend_title="Samples")


    with open(output_file, 'w', encoding='utf-8') as g:
        g.write(figure.to_html(full_html=False, include_plotlyjs='cdn'))


coverage_plot(bam_ref_files=snakemake.input, output_file=snakemake.output[0])
