import pandas as pd
import sys

# logging
sys.stderr = open(snakemake.log[0], "w")

df = pd.read_csv(snakemake.input[0])
df.columns.values[2] = snakemake.params.samples
df = df[df.iloc[:, 0] == snakemake.params.state].iloc[:, [1, 2]].T
df.to_csv(snakemake.output[0], header=False)