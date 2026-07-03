bits=("b2048" "b1024" "b512" "b256")

for bit in "${bits[@]}"
do
    h5file="chembl_35_${bit}.h5"
    echo "Running benchmarks for ${bit}"
    python benchmarking.py ${h5file} > benchmark_${bit}_0_7.log
done