#!/bin/bash
set -e

# Make sure executable is built
make

# Verify binary was built
if [ ! -f "ampt" ]; then
    echo "ERROR: ampt binary not found after make. Aborting."
    exit 1
fi

# Array of configurations
names=("default" "modified" "density2" "density3")
rhos=("0.0" "1.0" "2.0" "3.0")
iqmcs=("0" "1" "1" "1")

# Clean old runs
for name in "${names[@]}"; do
    rm -rf "run_${name}"
done

mkdir -p ana

for i in "${!names[@]}"; do
    name="${names[$i]}"
    rho="${rhos[$i]}"
    iqmc="${iqmcs[$i]}"
    
    echo "Starting configuration: $name (rho: $rho, iqmc: $iqmc)"
    
    # Create isolated directory
    mkdir "run_${name}" && cp ampt input.ampt model_data.csv "run_${name}/"
    
    # Write input.density
    echo "$rho" > "run_${name}/input.density"
    echo "$iqmc" >> "run_${name}/input.density"
    
    # Execute in background
    (
        cd "run_${name}"
        mkdir ana && (date '+%d%H%M%S' > nseed_runtime) && ./ampt < nseed_runtime > "ampt_stdout.log" 2>&1
        if [ -f ana/ampt.dat ]; then
            cp ana/ampt.dat "../ana/ampt_${name}.dat"
            echo "Finished $name: successfully copied ampt.dat"
        else
            echo "Failed $name: ana/ampt.dat not found"
        fi
    ) &
done

echo "Waiting for all AMPT jobs to finish..."
wait
echo "All AMPT simulations completed successfully!"

# Kept for debugging
# for name in "${names[@]}"; do
#     rm -rf "run_${name}"
# done

echo "Running Python plotting scripts..."
python3 plot_pt_spectra.py
python3 plot_ratios.py
python3 plot_v1_v2.py
python3 plot_extra_observables.py
python3 plot_advanced_baryon_stopping.py
python3 plot_advanced_splittings.py

echo "Data production and plotting complete!"
