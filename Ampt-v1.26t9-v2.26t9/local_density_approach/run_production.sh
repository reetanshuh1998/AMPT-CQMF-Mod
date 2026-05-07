#!/bin/bash
# =============================================================
# Full Production Run: AMPT-CQMF Local Density Phase 1
# Au+Au @ 7.7 GeV, Minimum Bias (b=0-14 fm)
# Configurations: Default, Fixed Density (1,2,3 rho0), Local Density
# =============================================================
set -e

NUM_EVENTS=2000

names=("default" "fixed_rho1" "fixed_rho2" "fixed_rho3" "local_density")
rhos=("0.0"      "1.0"        "2.0"        "3.0"        "0.0")
iqmcs=("0"       "1"          "1"          "1"          "2")

mkdir -p ana

echo "============================================================"
echo " AMPT-CQMF Production: Local Density Phase 1"
echo " Energy: 7.7 GeV Au+Au | Events: $NUM_EVENTS | BMAX: 14 fm"
echo "============================================================"

for i in "${!names[@]}"; do
    name="${names[$i]}"
    rho="${rhos[$i]}"
    iqmc="${iqmcs[$i]}"

    rm -rf "run_${name}"
    mkdir "run_${name}"
    cp ampt input.ampt model_data.csv "run_${name}/"

    # Set number of events
    sed -i "s/^[0-9]*[[:space:]]*! NEVNT.*$/${NUM_EVENTS}    ! NEVNT (total number of events)/" "run_${name}/input.ampt"

    # Write density config
    echo "$rho"  > "run_${name}/input.density"
    echo "$iqmc" >> "run_${name}/input.density"

    echo "Launching: $name (rho=${rho} rho0, iqmc=${iqmc})"

    (
        cd "run_${name}"
        mkdir -p ana
        attempts=0
        success=0
        while [ $success -eq 0 ]; do
            attempts=$((attempts+1))
            SEED=$(($(date +%s) + $RANDOM + $i * 1000))
            echo "$SEED" > nseed_runtime
            ./ampt < nseed_runtime > ampt_stdout.log 2>&1
            if [ -f ana/ampt.dat ] && [ -s ana/ampt.dat ]; then
                cp ana/ampt.dat "../ana/ampt_${name}.dat"
                echo "SUCCESS: ${name} finished in ${attempts} attempt(s)"
                success=1
            else
                echo "RETRY [${name}]: attempt $attempts failed, retrying..."
                rm -rf ana && mkdir -p ana
            fi
        done
    ) &

done

echo ""
echo "All ${#names[@]} configurations launched in parallel."
echo "Monitoring: tail -f run_<name>/ampt_stdout.log"
echo ""

wait
echo "============================================================"
echo " ALL PRODUCTION RUNS COMPLETE"
echo "============================================================"
ls -lh ana/ampt_*.dat 2>/dev/null || echo "WARNING: Some .dat files missing!"
