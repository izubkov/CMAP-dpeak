#!/bin/bash

#################################################################################
## Sandbox initialization.
#
## Pulls a sample solution (current DPeak algo), and tags it as cmap/solution
#docker pull cmap/sig_dpeak_tool:v58daea2b626a5eccb6523f5019f1bbf3b9966d13
#docker tag cmap/sig_dpeak_tool:v58daea2b626a5eccb6523f5019f1bbf3b9966d13 cmap/solution
#
## Builds scorer container.
#docker build -t cmap/scorer ./scorer

################################################################################
# Execution of tests. Each of local testing, provisional testing, and final
# testing, consist of two test cases.

TEST_CASE_1="DPK.CP001_A549_24H_X1_B42"
TEST_CASE_2="LITMUS.KD017_A549_96H_X1_B42"

exec_test() {
  docker run --rm -it \
    -v $(pwd)/input:/input \
    -v $(pwd)/output:/output \
    cmap/solution \
      --dspath /input/$1 \
      --out /output \
      --create_subdir 0 \
      --plate $1
}

exec_R() {
  Rscript exec.R \
    $1 \
    output
}

# Cleans output folder.
#
# TODO: docker scorer somehow writes as root
#
sudo rm ./output/* -rf

START_TIME=`date +%s`
#exec_test $TEST_CASE_1
exec_R $TEST_CASE_1
MID_TIME=`date +%s`
#exec_test $TEST_CASE_2
exec_R $TEST_CASE_2
END_TIME=`date +%s`

TEST_TIME_1=$((MID_TIME-START_TIME))
TEST_TIME_2=$((END_TIME-MID_TIME))

echo "Test #1 time [sec] =" $TEST_TIME_1
echo "Test #2 time [sec] =" $TEST_TIME_2

################################################################################
# Scoring.

docker run --rm \
  -v $(pwd)/output:/workdir \
  -v $(pwd)/ground-truth:/ground-truth \
  cmap/scorer $TEST_CASE_1 $TEST_CASE_2 $TEST_TIME_1 $TEST_TIME_2
