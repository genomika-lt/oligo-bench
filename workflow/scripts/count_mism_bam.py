import sys
import pysam
import pandas as pd

# logging
sys.stderr = open(snakemake.log[0], "w")

def count_mism_bam(bam, out, threads):
    
    pysam.index(bam)
    
    inp = pysam.AlignmentFile(bam, 'rb', threads=threads)
    
    # init a dict to count mismaches
    mism_counts = {
        'NM0': 0, 'NM1': 0, 'NM2': 0,
        'NM3': 0, 'NM4': 0, '>NM5': 0
    }
    
    # parse the BAM file
    for read in inp.fetch():
        # use primary reads
        if not read.is_supplementary and not read.is_secondary:
            nm = read.get_tag("NM")
            if nm == 0:
                mism_counts['NM0'] += 1
            elif nm == 1:
                mism_counts['NM1'] += 1
            elif nm == 2:
                mism_counts['NM2'] += 1
            elif nm == 3:
                mism_counts['NM3'] += 1
            elif nm == 4:
                mism_counts['NM4'] += 1
            else:
                mism_counts['>NM5'] += 1
    
    inp.close()
    
    # flatten
    mism_df = pd.DataFrame(
        list(mism_counts.items()), columns=['mism', 'count']
    )

    # write output
    mism_df.to_csv(out, header=False, index=False)


count_mism_bam(
    bam = snakemake.input[0],
    out = snakemake.output[0],
    threads = snakemake.threads
)
