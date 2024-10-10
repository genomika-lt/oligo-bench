#!/bin/bash
echo "SNAKEMAKE"
snakemake --lint
printf "\n\n"

echo "PYLINT"
pylint $(git ls-files *.py)
printf "\n\n"

echo "SNAKEFMT"
snakefmt ./