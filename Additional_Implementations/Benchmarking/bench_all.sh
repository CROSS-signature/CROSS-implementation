#!/bin/bash

for i in `ls build/bin/CROSS_benchmark*`
do
  # CPU 0 is a P-core on Intel i7-12700
  echo Benchmarking $i
  taskset --cpu-list 0 $i -T 2>&1 |grep TIME >> log.txt
done
