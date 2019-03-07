#!/bin/bash

date=$(date "+%Y_%m_%d")
file="snakemake_logs/${date}_multiblup.log"

snakemake --jobs 1000 --rerun-incomplete --latency-wait 60 --cluster-config submit.json --cluster "sbatch -p {cluster.p} -o {cluster.o} --cpus-per-task {cluster.cpus-per-task}" -p &>> ${file}
