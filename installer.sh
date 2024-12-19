#!/bin/bash
sudo echo "Started Installation"

wget https://github.com/genomika-lt/oligo-bench/archive/refs/heads/main.zip
unzip main.zip
rm main.zip

mv ./oligo-bench-main ./oligo_bench
sudo chmod -R a+rwx ./oligo_bench

wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
sudo chmod a+rwx Miniforge3-Linux-x86_64.sh
./Miniforge3-Linux-x86_64.sh -b
rm ./Miniforge3-Linux-x86_64.sh
~/miniforge3/condabin/conda init --all

echo "Installing packages"
~/miniforge3/condabin/conda create -c conda-forge -c bioconda -p oligo_bench/oligo snakemake minimap2 last samtools pandas snakefmt pysam plotly PyYAML requests -y
echo "Downloading dorado"
wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.9.0-linux-x64.tar.gz
echo "Unzipping dorado archive"
tar -xvzf dorado-0.9.0-linux-x64.tar.gz
rm -rf dorado-0.9.0-linux-x64.tar.gz
mv ./dorado-0.9.0-linux-x64 ./dorado
sudo chmod a+rwx ./dorado/bin/dorado
