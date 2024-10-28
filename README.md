# Snakemake workflow: `GENQC`

[![Snakemake](https://img.shields.io/badge/snakemake->=6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/jsimonas/oligo-bench/workflows/Tests/badge.svg)](https://github.com/jsimonas/oligo-bench/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake-based tool for quality control of synthetic oligonucleotide sequencing data produced by ONT sequencers.

## Installation

Snakemake is required to run this workflow. Installation of Snakemake can be found <a href='https://snakemake.readthedocs.io/en/stable/getting_started/installation.html'>here</a>

# Usage

```bash
snakemake --cores 1 --use-conda
```