import pandas as pd
import sys

# logging
sys.stderr = open(snakemake.log[0], "w")


def parse_samtools_stats(inp, samples, out):
    
    mq_lines = []
    id_lines = []
    co_lines = []
    
    with open(inp, 'r') as f:
        for line in f:
            if line.startswith("MAPQ"):
                mq_lines.append(line)
            elif line.startswith("ID"):
                id_lines.append(line)
            elif line.startswith("COV"):
                co_lines.append(line)

    mq_df = pd.DataFrame(
        [x.strip().split('\t') for x in mq_lines],
        columns=["Type", "MAPQ", samples]
    ).iloc[:, [2]].T
    
    in_df = pd.DataFrame(
        [x.strip().split('\t') for x in id_lines],
        columns=["Type", "Length", samples, "Num_Del"]
    ).iloc[:, [2]].T
    
    dl_df = pd.DataFrame(
        [x.strip().split('\t') for x in id_lines],
        columns=["Type", "Length", "Num_Ins", samples]
    ).iloc[:, [3]].T
        
    co_df = pd.DataFrame(
        [x.strip().split('\t') for x in co_lines],
        columns=["Type", "Range", "Depth", samples]
    ).iloc[:, [3]].T

    mq_df.to_csv(
        out[0],
        header=False
    )
    in_df.to_csv(
        out[1],
        header=False
    )
    dl_df.to_csv(
        out[2],
        header=False
    )
    co_df.to_csv(
        out[3],
        header=False
    )


parse_samtools_stats(
    inp = snakemake.input[0],
    samples = snakemake.params.samples,
    out = snakemake.output
)
