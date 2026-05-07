#!/bin/bash
# run_ampt_chunk.sh
# Usage: ./run_ampt_chunk.sh <iqmc_mode> <config_name> <rho>
#   iqmc_mode: 0 = Default, 2 = Local Density

IQMC=$1
CONFIG_NAME=$2
RHO=$3
if [ -z "$RHO" ]; then
    echo "Usage: $0 <iqmc_mode> <config_name> <rho>"
    exit 1
fi

EVENTS_PER_JOB=2000
JOB_INDEX=${LSB_JOBINDEX:-1}

# Define base directory (assuming script is run from inside kekcc_submission)
BASE_DIR=$(pwd)
OUTPUT_DIR="${BASE_DIR}/ana/chunks"
mkdir -p "${OUTPUT_DIR}"

# Create a unique working directory under group storage (avoid /tmp quota issues)
TMP_DIR="${BASE_DIR}/tmp_${CONFIG_NAME}_${JOB_INDEX}_$$"
mkdir -p "${TMP_DIR}"
trap "rm -rf ${TMP_DIR}" EXIT
cd "${TMP_DIR}"

# Copy necessary files
cp "${BASE_DIR}/ampt" .
cp "${BASE_DIR}/model_data.csv" .
cp "${BASE_DIR}/input.ampt.template" ./input.ampt
cp "${BASE_DIR}/ampt_to_root.py" .

# Modify input.ampt for this chunk (100 events, unique seed setup)
sed -i "s/^[0-9]\+[[:space:]]*![[:space:]]*NEVNT/${EVENTS_PER_JOB}    ! NEVNT/g" input.ampt

# Setup input.density
echo "${RHO}" > input.density
echo "${IQMC}" >> input.density

# Generate a unique random seed for this job that fits in a 32-bit integer
SHORT_JOBID=$(echo ${LSB_JOBID:-1000} | tail -c 6)
UNIQUE_SEED=$(( ${SHORT_JOBID} * 10000 + ${JOB_INDEX} + ${IQMC} * 12345 ))
echo "${UNIQUE_SEED}" > nseed_runtime

mkdir -p ana

# Run AMPT
./ampt < nseed_runtime > ampt.log 2>&1

# Convert to ROOT and immediately delete raw data to free /tmp space
python3 ampt_to_root.py ana/ampt.dat ana/ampt.root 2>&1
rm -f ana/ampt.dat

# Copy the generated data back to the central storage
if [ -s ana/ampt.root ]; then
    CHUNK_STR=$(printf "%03d" ${JOB_INDEX})
    cp ana/ampt.root "${OUTPUT_DIR}/ampt_${CONFIG_NAME}_chunk_${CHUNK_STR}.root"
    echo "Successfully generated chunk ${JOB_INDEX} for ${CONFIG_NAME}."
else
    echo "Error: ampt.root not found or empty."
    cat ampt.log 2>/dev/null
    exit 1
fi

# Cleanup
rm -rf "${TMP_DIR}"
