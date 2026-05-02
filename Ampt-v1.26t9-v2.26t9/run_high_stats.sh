#!/bin/bash
set -e

# Production Configuration
NUM_EVENTS=20000

# Array of configurations
names=("default" "modified" "density2" "density3")
rhos=("0.0" "1.0" "2.0" "3.0")
iqmcs=("0" "1" "1" "1")

# Clean old runs (optional, but safer to start fresh)
for name in "${names[@]}"; do
    rm -rf "run_${name}_highstats"
done

mkdir -p ana_highstats

echo "Starting High Statistics AMPT Production (20k events per config)"
echo "Target Energy: 7.7 GeV Au+Au"

for i in "${!names[@]}"; do
    name="${names[$i]}"
    rho="${rhos[$i]}"
    iqmc="${iqmcs[$i]}"
    
    echo "Launching: $name (rho: $rho, iqmc: $iqmc)"
    
    # Create isolated directory
    mkdir "run_${name}_highstats" && cp ampt input.ampt model_data.csv "run_${name}_highstats/"
    
    # Set the requested number of events dynamically
    sed -i "s/^[0-9]*[ \t]*! NEVNT.*$/${NUM_EVENTS}    ! NEVNT (total number of events)/" "run_${name}_highstats/input.ampt"
    
    # Write input.density
    echo "$rho" > "run_${name}_highstats/input.density"
    echo "$iqmc" >> "run_${name}_highstats/input.density"
    
    # Execute in background
    (
        cd "run_${name}_highstats"
        mkdir -p ana 
        # Generate seed based on time + index to avoid collisions
        SEED=$(($(date +%s) + $i))
        echo "$SEED" > nseed_runtime
        ./ampt < nseed_runtime > "ampt_stdout.log" 2>&1
        
        if [ -f ana/ampt.dat ]; then
            cp ana/ampt.dat "../ana_highstats/ampt_${name}.dat"
            echo "SUCCESS: $name production finished and copied to ana_highstats/"
        else
            echo "ERROR: $name production failed (ana/ampt.dat not found)"
        fi
    ) &
done

echo "Waiting for all high-statistics jobs to finish (approx 100-120 minutes)..."
wait
echo "All high-statistics AMPT simulations completed!"
