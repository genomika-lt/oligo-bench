
rule count_total_reads_number:
    input:
        expand(
            "results/basecalled/{sample_id}.bam",
            sample_id=samples["sample_id"],
        ),
    output:
        "results/statistics/total/count_total_reads_number.csv",
    log:
        "logs/statistics/overall/count_total_reads_number.log",
    conda:
        "../../envs/pysam.yaml"
    script:
        "../../scripts/statistics/overall/summary_table/count_total_reads_number.py"


rule count_total_bases_number:
    input:
        expand(
            "results/basecalled/{sample_id}.bam",
            sample_id=samples["sample_id"],
        ),
    output:
        "results/statistics/total/count_total_bases_number.csv",
    log:
        "logs/statistics/overall/count_total_bases_number.log",
    conda:
        "../../envs/pysam.yaml"
    script:
        "../../scripts/statistics/overall/summary_table/count_total_bases_number.py"


rule count_passed_reads_number:
    input:
        expand(
            "results/basecalled/passed_{sample_id}.bam",
            sample_id=samples["sample_id"],
        ),
    output:
        "results/statistics/total/count_passed_reads_number.csv",
    log:
        "logs/statistics/overall/count_passed_reads_number.log",
    conda:
        "../../envs/pysam.yaml"
    script:
        "../../scripts/statistics/overall/summary_table/count_passed_reads_number.py"


rule count_aligned_reads_number:
    input:
        expand(
            "results/aligned/passed_{sample_id}.bam",
            sample_id=samples["sample_id"],
        ),
    output:
        "results/statistics/total/count_aligned_reads_number.csv",
    log:
        "logs/statistics/overall/count_aligned_reads_number.log",
    conda:
        "../../envs/pysam.yaml"
    script:
        "../../scripts/statistics/overall/summary_table/count_aligned_reads_number.py"


rule count_passed_bases_number:
    input:
        expand(
            "results/basecalled/passed_{sample_id}.bam",
            sample_id=samples["sample_id"],
        ),
    output:
        "results/statistics/total/count_passed_bases_number.csv",
    log:
        "logs/statistics/overall/count_passed_bases_number.log",
    conda:
        "../../envs/pysam.yaml"
    script:
        "../../scripts/statistics/overall/summary_table/count_passed_bases_number.py"


rule count_gc_passed_bases_number:
    input:
        expand(
            "results/basecalled/passed_{sample_id}.bam",
            sample_id=samples["sample_id"],
        ),
    output:
        "results/statistics/total/count_gc_passed_bases_number.csv",
    log:
        "logs/statistics/overall/count_gc_passed_bases_number.log",
    conda:
        "../../envs/pysam.yaml"
    script:
        "../../scripts/statistics/overall/summary_table/count_gc_passed_bases_number.py"


rule experiment_timestamps:
    input:
        samples["path_to_sample"],
    output:
        "results/statistics/total/experiment_timestamps.csv",
    log:
        "logs/statistics/overall/experiment_timestamps.log",
    conda:
        "../../envs/pysam.yaml"
    script:
        "../../scripts/statistics/overall/summary_table/experiment_timestamps.py"


rule summary_table:
    input:
        total_reads="results/statistics/total/count_total_reads_number.csv",
        total_bases="results/statistics/total/count_total_bases_number.csv",
        passed_reads="results/statistics/total/count_passed_reads_number.csv",
        passed_bases="results/statistics/total/count_passed_bases_number.csv",
        gc_passed_bases="results/statistics/total/count_gc_passed_bases_number.csv",
        timestamps="results/statistics/total/experiment_timestamps.csv",
    output:
        "results/statistics/summary_table.html",
    log:
        "logs/summary_table.log",
    conda:
        "../../envs/plotly.yaml"
    script:
        "../../scripts/statistics/overall/summary_table/summary_table.py"


rule filter_2_primers:
    input:
        "results/mapping/{sample_id}.sam",
    output:
        "results/mapping/2_primers_reads_{sample_id}.sam",
    log:
        "logs/mapping/2_primers_reads_{sample_id}.log",
    conda:
        "../../envs/pysam.yaml"
    script:
        "../../scripts/statistics/overall/mapping_table/filter_2_primers.py"


rule read_quality_histogram:
    input:
        expand(
            "results/basecalled/{sample_id}.bam",
            sample_id=samples["sample_id"],
        ),
    output:
        "results/statistics/read_quality_histogram.html",
    log:
        "logs/read_quality_histogram.log",
    conda:
        "../../envs/plotly.yaml"
    script:
        "../../scripts/statistics/overall/read_quality_histogram.py"


rule passed_read_length_histogram:
    input:
        expand(
            "results/basecalled/passed_{sample_id}.bam",
            sample_id=samples["sample_id"],
        ),
    output:
        "results/statistics/passed_read_length_histogram.html",
    log:
        "logs/passed_read_length_histogram.log",
    conda:
        "../../envs/plotly.yaml"
    script:
        "../../scripts/statistics/overall/passed_read_length_histogram.py"


rule quality_per_base_position:
    input:
        expand(
            "results/basecalled/passed_{sample_id}.bam",
            sample_id=samples["sample_id"],
        ),
    output:
        "results/statistics/quality_per_base_position.html",
    log:
        "logs/quality_per_base_position.log",
    conda:
        "../../envs/plotly.yaml"
    script:
        "../../scripts/statistics/overall/quality_per_base_position.py"


rule primers_number_histogram:
    input:
        expand("results/mapping/{sample_id}.sam", sample_id=samples["sample_id"]),
    output:
        "results/statistics/primers_number_histogram.html",
    log:
        "logs/primers_number_histogram.log",
    conda:
        "../../envs/pysam.yaml"
    script:
        "../../scripts/statistics/overall/primers_number_histogram.py"


rule count_mapped_primers_reads_number:
    input:
        expand("results/mapping/2_primers_reads_{sample_id}.sam", sample_id=samples["sample_id"]),
    output:
        "results/statistics/total/mapped_primers_reads_number.csv",
    log:
        "logs/count_mapped_primers_reads_number.log",
    conda:
        "../../envs/pysam.yaml"
    script:
        "../../scripts/statistics/overall/mapping_table/count_mapped_primers_reads_number.py"


rule count_mapped_primers_bases_number:
    input:
        expand("results/mapping/2_primers_reads_{sample_id}.sam", sample_id=samples["sample_id"]),
    output:
        "results/statistics/total/mapped_primers_bases_number.csv",
    log:
        "logs/count_mapped_primers_bases_number.log",
    conda:
        "../../envs/pysam.yaml"
    script:
        "../../scripts/statistics/overall/mapping_table/count_mapped_primers_bases_number.py"


rule count_forward_mapped_reads_percentage:
    input:
        expand("results/mapping/2_primers_reads_{sample_id}.sam", sample_id=samples["sample_id"]),
    output:
        "results/statistics/total/forward_mapped_reads_number.csv",
    log:
        "logs/mapping/count_forward_mapped_reads_percentage.log"
    conda:
        "../../envs/pysam.yaml"
    script:
        "../../scripts/statistics/overall/mapping_table/count_forward_mapped_reads_percentage.py"


rule mapping_table:
    input:
        forward_reads="results/statistics/total/forward_mapped_reads_number.csv",
        mapped_reads="results/statistics/total/mapped_primers_reads_number.csv",
        mapped_bases="results/statistics/total/mapped_primers_bases_number.csv",
        total_bases="results/statistics/total/count_total_bases_number.csv",
        total_reads="results/statistics/total/count_total_reads_number.csv",
        aligned_reads="results/statistics/total/count_aligned_reads_number.csv",
    output:
        "results/statistics/mapping_table.html",
    log:
        "logs/mapping_table.log",
    conda:
        "../../envs/plotly.yaml"
    script:
        "../../scripts/statistics/overall/mapping_table/mapping_table.py"


rule payload_errors_number:
    input:
        expand("results/aligned/passed_{sample_id}.bam", sample_id=samples["sample_id"]),
        expand("results/mapping/primers_{sample_id}.fa", sample_id=samples["sample_id"]),
    output:
        "results/statistics/total/payload_errors_number.csv"
    log:
        "logs/payload_errors_number.log",
    conda:
        "../../envs/plotly.yaml"
    script:
        "../../scripts/statistics/overall/payload_table/payload_errors_number.py"



rule payload_table:
    input:
        payload_errors="results/statistics/total/payload_errors_number.csv",
        total_reads="results/statistics/total/count_total_reads_number.csv",
    output:
        "results/statistics/payload_table.html",
    log:
        "logs/payload_table.log",
    conda:
        "../../envs/plotly.yaml"
    script:
        "../../scripts/statistics/overall/payload_table/payload_table.py"


rule coverage_plot:
    input:
        expand(
            "results/aligned/passed_{sample_id}.bam",
            sample_id=samples["sample_id"],
        ),
        expand(
            "results/reference/{ref_name}.fa",
            ref_name=samples["ref_name"],
        ),
    output:
        "results/statistics/coverage_plot.html",
    log:
        "logs/statistics/coverage_plot.log"
    conda:
        "../../envs/plotly.yaml"
    script:
        "../../scripts/statistics/overall/coverage_plot.py"
