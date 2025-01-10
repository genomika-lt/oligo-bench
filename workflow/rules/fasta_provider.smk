rule provide_fasta:
    params:
         input_txt = get_reference_name
    output:
        "results/reference/{ref_name}.fa"
    log:
        "logs/convert_txt_fasta_{ref_name}.log"
    script:
        "../scripts/txt_fasta_convert.py"

