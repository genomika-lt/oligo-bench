import pandas as pd
import sys

# logging
sys.stderr = open(snakemake.log[0], "w")

def parse_pore_activity(inp, samples, state, out):
    df = pd.read_csv(inp)
    df.columns.values[2] = samples
    df = df[df.iloc[:, 0] == state].iloc[:, [2]].T
    df.to_csv(out, header=False)


parse_pore_activity(
    inp = snakemake.input[0],
    samples = snakemake.params.samples,
    state = snakemake.params.state,
    out = snakemake.output[0]
)
