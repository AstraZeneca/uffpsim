bits=("256" "512" "1024" "2048")
CWD=$(pwd)

# First create the initial database with cutoff 0.2 for all sizes
for size in "${bits[@]}"
do
    echo "Creating database for size $size"
    python create_db.py $size
done

# Re-do clustering for each size with the specified cutoffs
for size in "${bits[@]}"
do
    echo "Re-doing clustering for size $size"
    cd "$CWD/b$size"
    python redo_clustering.py $size
done
cd "$CWD"

