import pandas as pd
import sys

# logging
sys.stderr = open(snakemake.log[0], "w")

def add_id_to_csv(inp, col_id, samples, out):
    
    df = pd.read_csv(inp)
    df.insert(0, col_id, samples)
    df.to_csv(out, index=False)


add_id_to_csv(
    inp = snakemake.input[0],
    col_id = snakemake.params.id,
    samples = snakemake.params.samples,
    out = snakemake.output[0]
)