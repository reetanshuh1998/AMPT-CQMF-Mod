#!/bin/bash
# merge_chunks.sh
# Combines all chunk files into the final datasets

BASE_DIR=$(pwd)
CHUNK_DIR="${BASE_DIR}/ana/chunks"
OUTPUT_DIR="${BASE_DIR}/ana"

mkdir -p "${OUTPUT_DIR}"

CONFIGS=("default" "localdensity" "fixed_rho1" "fixed_rho2" "fixed_rho3")

for CONF in "${CONFIGS[@]}"; do
    if ls "${CHUNK_DIR}"/ampt_${CONF}_chunk_*.root 1> /dev/null 2>&1; then
        echo "Merging ${CONF} AMPT chunks..."
        rm -f "${OUTPUT_DIR}/ampt_${CONF}.root"
        hadd -f -O "${OUTPUT_DIR}/ampt_${CONF}.root" "${CHUNK_DIR}"/ampt_${CONF}_chunk_*.root
        echo "  -> Created ana/ampt_${CONF}.root"
    fi
done

echo "Done!"
