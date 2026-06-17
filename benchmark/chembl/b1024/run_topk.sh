cutoff=(0.7 0.6 0.5)
topk=(1 10 100 1000)
for c in "${cutoff[@]}"
do
    for k in "${topk[@]}"
    do
        echo "Running benchmark for search-cutoff ${c} and topk ${k}"
        python benchmarking_topk.py chembl35_18.h5 $c $k > topk/benchmark_${c}_${k}.log
    done
done
