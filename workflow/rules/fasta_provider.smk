rule provide_fasta:
    params:
         input_txt = get_reference
    output:
        "results/reference/{sample_id}.fa"
    log:
        "logs/convert_txt_fasta_{sample_id}.log"
    script:
        "../scripts/txt_fasta_convert.py"

