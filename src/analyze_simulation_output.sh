#!/bin/bash
FILES=../data/data_dump/forbeta/*
let COUNTER=0
for f in $FILES;
do
  let COUNTER=COUNTER+1
  for i in `seq 1 4`;
  do
    ../stan/reflected_beta sample num_samples=4000 num_warmup=4000 thin=1 \
      id=$i \
      random seed=420 \
      data file=$f \
      output file=../data/mcmc_out/refbeta_wei_${COUNTER}_${i}.csv &
  done
  wait
  for i in `seq 1 4`;
  do
    ../stan/refbeta_exp sample num_samples=4000 num_warmup=4000 thin=1 \
      id=$i \
      random seed=420 \
      data file=$f \
      output file=../data/mcmc_out/refbeta_exp_${COUNTER}_${i}.csv &
  done
  wait
  for i in `seq 1 4`;
  do
    ../stan/weibull sample num_samples=4000 num_warmup=4000 thin=1 \
      id=$i \
      random seed=420 \
      data file=$f \
      output file=../data/mcmc_out/weibull_${COUNTER}_${i}.csv &
  done
  wait
  for i in `seq 1 4`;
  do
    ../stan/exponential sample num_samples=4000 num_warmup=4000 thin=1 \
      id=$i \
      random seed=420 \
      data file=$f \
      output file=../data/mcmc_out/exponential_${COUNTER}_${i}.csv &
  done
  wait
done

FILES=../data/data_dump/forpert/*
NEWCOUNTER=0
for f in $FILES;
do
  let NEWCOUNTER=NEWCOUNTER+1
  for i in `seq 1 4`;
  do
    ../stan/pert sample num_samples=4000 num_warmup=4000 thin=1 \
      id=$i \
      random seed=420 \
      data file=$f \
      output file=../data/mcmc_out/pert_wei_${NEWCOUNTER}_${i}.csv &
  done
  wait
  for i in `seq 1 4`;
  do
    ../stan/pert_exp sample num_samples=4000 num_warmup=4000 thin=1 \
      id=$i \
      random seed=420 \
      data file=$f \
      output file=../data/mcmc_out/pert_exp_${NEWCOUNTER}_${i}.csv &
  done
  wait
done
