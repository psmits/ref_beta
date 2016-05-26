#!/bin/bash
FILES=../data/data_dump/special*
for f in $FILES;
do
  for i in `seq 1 4`;
  do
    ../stan/test/reflected_beta sample num_samples=1000 num_warmup=1000 thin=1 \
      id=$i \
      random seed=420 \
      data file=$f \
      output file=../data/mcmc_out/pert_attempts_${i}.csv &
  done
  wait
done
