thresholds=(10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30)
for threshold in "${thresholds[@]}"
do
    echo "Running benchmark for clustering threshold ${threshold}"
    python benchmarking.py chembl35_${threshold}.h5 > benchmark_${threshold}.log
done

