rule minimap2:
    input:
        target=get_reference,
        query="results/basecalled/passed_{sample_id}.bam",
    output:
        "results/aligned/{sample_id}.bam",
    log:
        "logs/minimap2_{sample_id}.log",
    params:
        extra="-ax sr -L --MD -Y -k10 -w5 -m10",
        sorting="coordinate",
        sort_extra="",
    threads: 4
    wrapper:
        "v3.13.3/bio/minimap2/aligner"
