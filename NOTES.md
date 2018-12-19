Baseline:

Test #1 time [sec] = 163
Test #2 time [sec] = 161

DPK:

Median Spearman correlation: 0.701
AUC (hi and lo): 0.92 and 0.921
Overall accuracy score 1e6 * COR * AUC = 645653.205832184

LIT:

Median Spearman correlation: 0.819
AUC (hi and lo): 0.904 and 0.926
Overall accuracy score 1e6 * COR * AUC = 749249.376560217

OVERALL SCORE = 421009.415034549



Running
 > nohup  /cmap/tools/sig_tools/run_sig_dpeak_tool.sh /cmap/tools/sig_tools/mcr/versions/v84 --dspath /input/DPK.CP001_A549_24H_X1_B42 --out /output --create_subdir 0 --plate DPK.CP001_A549_24H_X1_B42 


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
