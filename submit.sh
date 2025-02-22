#!/bin/bash

date=$(date "+%Y_%m_%d")
file="snakemake_logs/${date}_multiblup.log"

snakemake -j 1000 \
--rerun-incomplete \
--latency-wait 60 \
--cluster-config submit.json \
--cluster "sbatch -p {cluster.p} -o {cluster.o} --job-name {cluster.name} --cpus-per-task {cluster.cpus-per-task}"
