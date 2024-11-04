rule extract_primers:
    params:
        forward_primer=get_forward_primer,
        reverse_primer=get_reverse_primer,
    output:
        "results/mapping/primers_{sample_id}.fa",
    log:
        "logs/extract_primers_{sample_id}.log",
    conda:
        "../../envs/pysam.yaml"
    script:
        "../../scripts/extract_primers.py"


rule last:
    input:
        reference_file="results/mapping/primers_{sample_id}.fa",
        reads="results/basecalled/passed_{sample_id}.fastq",
    output:
        aligned="results/mapping/{sample_id}.maf",
        converted="results/mapping/{sample_id}.sam",
    log:
        "logs/last_aligner_{sample_id}.log",
    params:
        database="results/mapping/{sample_id}",
        train="results/mapping/{sample_id}.train",
    threads: 4
    conda:
        "../../envs/pysam.yaml"
    shell:
        """
        lastdb -P {threads} {params.database} {input.reference_file} && \
        last-train -s 2 -P {threads} {params.database} -Q 1 {input.reads} > {params.train} && \
        lastal -s 2 -P {threads} -Q 1 -p {params.train} {params.database} {input.reads} > {output.aligned} && \
        maf-convert sam -d {output.aligned} > {output.converted}
        """
