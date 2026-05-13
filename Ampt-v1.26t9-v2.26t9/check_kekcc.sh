echo "Checking chunks:"
ls -lh kekcc_20x20x20_production/kekcc_submission/ana/chunks | wc -l
echo "Checking logs for errors:"
grep -i "error" kekcc_20x20x20_production/kekcc_submission/logs/err_*.log | head -n 5
