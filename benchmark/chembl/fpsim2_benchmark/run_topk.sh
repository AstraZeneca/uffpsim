bits=("b2048" "b1024" "b512" "b256")
K=(1 10 100 1000)

for bit in "${bits[@]}"
do
    h5file="chembl_35_${bit}.h5"
    echo "Running benchmarks for ${bit}"

    for k in "${K[@]}"
    do
        echo "t=0.7 topk=$k"
        python benchmarking_topk.py ${h5file} 0.7 $k > topk/benchmark_${bit}_7_${k}.log

        echo "t=0.6 topk=$k"
        python benchmarking_topk.py ${h5file} 0.6 $k > topk/benchmark_${bit}_6_${k}.log

        echo "t=0.5 topk=$k"
        python benchmarking_topk.py ${h5file} 0.5 $k > topk/benchmark_${bit}_5_${k}.log
    done
done
