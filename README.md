# Snakemake workflow: `GENQC`

[![Snakemake](https://img.shields.io/badge/snakemake->=6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/jsimonas/oligo-bench/workflows/Tests/badge.svg)](https://github.com/jsimonas/oligo-bench/actions?query=branch%3Amain+workflow%3ATests)
[![PyQt6](https://img.shields.io/badge/PyQt6-%3E%3D6.0-brightgreen.svg)](https://pypi.org/project/PyQt6/)

A Snakemake-based tool for quality control of synthetic oligonucleotide sequencing data produced by ONT sequencers.

## Installation

Snakemake is required to run this workflow. Installation of Snakemake can be found [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). 

In addition to Snakemake, this project includes a PyQt6 graphical user interface (GUI) application for generating YAML and CSV files needed for analyses. 

Ensure you have Python installed. You can download it from the [official Python website](https://www.python.org/downloads/) and remember to check "Add Python to PATH" during installation.
PyQt6 can be installed via pip:
```bash
pip3 install PyQt6
```

## Usage

### Running the Snakemake Workflow
To run the Snakemake workflow, use the following command:

```bash
snakemake --cores 1 --use-conda
```
### Launching the GUI Application
The GUI application, implemented in PyQt6, provides an interface for creating YAML and CSV files.

To run the GUI application navigate to the folder `config` containing `config_gui.bin` and execute it.
