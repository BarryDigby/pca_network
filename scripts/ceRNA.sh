#!/usr/bin/env bash

set -e

# run rscripts to generate ceRNA network results
# subsequent R plots are for downstream analysis of ceRNA net

for i in `seq 1 1 7`; do
  echo "starting script ${i}"
  Rscript ${i}.*.R
done
