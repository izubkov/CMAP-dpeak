#!/bin/bash

cd submission

echo "Starting exec.R ..."
echo $1
echo $2
echo $(pwd)
Rscript exec.R $1 $2 $3 $4 $5 $6 $7 $8 $9
