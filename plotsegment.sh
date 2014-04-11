for i in {0..8}
do
    echo "plotting $i"
    doplot $1 $((i*1000000)) 1000000
done
