#!/bin/bash
sudo echo "Started Installation"

wget https://github.com/genomika-lt/oligo-bench/archive/refs/heads/main.zip && \
unzip main.zip && \
rm main.zip && \

mv ./oligo-bench-main ./oligo_bench && \
sudo chmod -R a+rwx ./oligo_bench

wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh && \
sudo chmod a+rwx Miniforge3-Linux-x86_64.sh && \
./Miniforge3-Linux-x86_64.sh -b && \
rm ./Miniforge3-Linux-x86_64.sh && \
~/miniforge3/condabin/conda init --all

echo "Installing packages"
~/miniforge3/condabin/conda create -c conda-forge -c bioconda -p ./oligo snakemake minimap2 last samtools pandas snakefmt pysam plotly PyYAML requests -y