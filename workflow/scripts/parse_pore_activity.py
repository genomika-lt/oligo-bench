import pandas as pd
import sys

# logging
sys.stderr = open(snakemake.log[0], "w")

def parse_pore_activity(inp, samples, out):
    
    df = pd.read_csv(inp)
    
    df.columns.values[2] = samples
    
    act_df = df[df.iloc[:, 0] == "pore"].iloc[:, [2]].reset_index(drop=True).T
    seq_df = df[df.iloc[:, 0] == "strand"].iloc[:, [2]].reset_index(drop=True).T
    
    act_df.replace(0, 1, inplace=True)
    seq_df.replace(0, 1, inplace=True)
    
    if act_df.shape != seq_df.shape:
        raise ValueError("pore_act and pore_seq have different shapes")
    
    ratio = (seq_df / act_df) * 100
        
    act_df.to_csv(out[0], header=False)
    seq_df.to_csv(out[1], header=False)
    ratio.to_csv(out[2], header=False)


parse_pore_activity(
    inp = snakemake.input[0],
    samples = snakemake.params.samples,
    out = snakemake.output
)
