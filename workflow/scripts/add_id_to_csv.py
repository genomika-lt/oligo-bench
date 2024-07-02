import pandas as pd
import sys

# logging
sys.stderr = open(snakemake.log[0], "w")

df = pd.read_csv(snakemake.input[0])
df.insert(0, snakemake.params.id, snakemake.params.samples)
df.to_csv(snakemake.output[0], index=False)