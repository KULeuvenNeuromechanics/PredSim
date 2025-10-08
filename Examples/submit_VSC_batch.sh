#!/bin/bash

# first make the file executable: $ chmod +x Examples/submit_VSC_batch.sh
# run this script: $ bash Examples/submit_VSC_batch.sh

subjects=("Falisse_et_al_2022")

for subject in "${subjects[@]}"; do
  for velocity in $(seq 1.0 0.1 1.2); do
    echo "Submitting job for $subject at $velocity m/s"
    sbatch --export=ALL,PREDSIM_SUBJECT=$subject,PREDSIM_VELOCITY=$velocity run_simulation_batch.slurm
  done
done