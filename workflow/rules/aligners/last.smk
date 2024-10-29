rule extract_primers:
    params:
        forward_primer=get_forward_primer,
        reverse_primer=get_reverse_primer
    output:
        config["output_folder"] + "/mapping/primers_{sample_id}.fa"
    log:
        "logs/extract_primers_{sample_id}.log"
    script:
        "../../scripts/extract_primers.py"


rule last:
    input:
        reference_file="primers_{sample_id}.fa",
        reads="passed_{sample_id}.fastq"
    output:
        "results/mapping/{sample_id}.bam",
    log:
        "logs/last_aligner_{sample_id}.log",
    shell:
        "lastdb {wildcards.sample_id} {input.reference_file}"
