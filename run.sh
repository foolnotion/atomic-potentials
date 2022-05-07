NAME=`echo $RANDOM | md5sum | head -c 20; echo;`

mkdir -p results/${NAME}/runs
cp "$(basename $0)" results/${NAME}/

for i in {1..50}
do
    SEED=`shuf -i 0-10000000000 -n 1`
    ./build/atomic-potentials \
        --dataset data/DFT_Cu/energy.data \
        --coordinates data/DFT_Cu/coord.data \
        --cluster-size 32 \
        --target energy \
        --reinserter keep-best \
        --evaluations 10000000 \
        --population-size 10000 \
        --pool-size 10000 \
        --female-selector tournament:17 \
        --epsilon 1e-5 \
        --iterations 0 \
        --train 0:75 \
        --test 75:150 \
        --enable-symbols pow \
        --maxlength 20 \
        --inputs r,q \
        --shuffle \
        --cutoff-radius 5.0 \
        --error-metric r2 \
        --generations 1000 \
        --seed ${SEED} \
        --results results/${NAME}/ | tee -a results/${NAME}/runs/${SEED}.csv
done
