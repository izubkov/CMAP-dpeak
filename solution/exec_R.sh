#!/bin/bash

echo "Starting exec.R ..."
echo $1
echo $2
echo $(pwd)
Rscript /solution/exec.R $2 $4
