#!/bin/bash
# submit_all.sh
# Master submission script for KEKCC LSF cluster

# Make scripts executable
chmod +x run_ampt_chunk.sh
chmod +x merge_chunks.sh

# Ensure logs and chunk output directories exist
mkdir -p logs
mkdir -p ana/chunks

# Number of chunks (100 chunks x 2000 events = 200,000 events total per config)
NUM_CHUNKS=100

# Define configurations: Name IQMC RHO
configs=(
    "default 0 0.0"
    "fixed_rho1 1 1.0"
    "fixed_rho2 1 2.0"
    "fixed_rho3 1 3.0"
    "localdensity 2 0.0"
)

for config in "${configs[@]}"; do
    read -r name iqmc rho <<< "$config"
    echo "Submitting ${name} (iqmc=${iqmc}, rho=${rho})..."
    bsub -q s -J "ampt_${name}[1-${NUM_CHUNKS}]%10" -o logs/out_${name}_%I.log -e logs/err_${name}_%I.log ./run_ampt_chunk.sh $iqmc $name $rho
done

echo "--------------------------------------------------------"
echo "Jobs submitted! You can monitor them using:"
echo "  bjobs"
echo ""
echo "Once ALL jobs are complete, merge the data by running:"
echo "  ./kekcc_submission/merge_chunks.sh"
echo "--------------------------------------------------------"
