#!/bin/bash
set -e
mkdir -p quick_ld && cd quick_ld
cp ../ampt ../model_data.csv ../input.ampt .
sed -i 's/^[0-9]*[ \t]*! NEVNT.*$/500    ! NEVNT (total number of events)/' input.ampt
echo "0.0" > input.density
echo "2" >> input.density
mkdir -p ana
success=0
attempts=0
while [ $success -eq 0 ]; do
    attempts=$((attempts+1))
    SEED=$(($(date +%s) + $RANDOM))
    echo "$SEED" > nseed_runtime
    ./ampt < nseed_runtime > ampt_stdout.log 2>&1
    if [ -f ana/ampt.dat ] && [ -s ana/ampt.dat ]; then
        cp ana/ampt.dat ../ana/ampt_localdensity_fixed.dat
        echo "SUCCESS: finished in ${attempts} attempt(s)"
        success=1
    else
        echo "RETRY: attempt $attempts failed, retrying..."
        rm -rf ana && mkdir -p ana
    fi
done
