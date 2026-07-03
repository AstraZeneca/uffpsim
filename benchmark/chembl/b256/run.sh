thresholds=(10 15 20 22 24 26 28 30 32 34 36 38 40 42 44)
for threshold in "${thresholds[@]}"
do
    echo "Running benchmark for clustering threshold ${threshold}"
    python benchmarking.py chembl35_${threshold}.h5 > benchmark_${threshold}.log
done

