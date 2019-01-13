## Publish new docker:

``` bash
cd hub-docker
export IMG_TAG=0
docker build -t izubkov/cmap:$IMG_TAG .
sudo docker push izubkov/cmap:$IMG_TAG
```

## Run locally:

``` bash
./run.sh
```

## Make submission:

``` bash
./run.sh
cd solution
# zip solution & submit
```

# Baseline:

```
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
```

### dpeak matlab routine:

constants and params (detect_lxb_peaks_single):

https://github.com/cmap/cmapM/blob/db40b0af5d3440c43bc4e8c6947399d9432cbc43/%2Bcmapm/%40Pipeline/private/detect_lxb_peaks_single.m#L49

ksdensity (kernel smoothing based on normal kernel function)

https://www.mathworks.com/help/stats/ksdensity.html

findpeaks

www.mathworks.com/help/signal/ref/findpeaks.html

UNI:

github.com/cmap/cmapM/blob/db40b0af5d3440c43bc4e8c6947399d9432cbc43/%2Bcmapm/%40Pipeline/private/dpeak_pipe.m#L62




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

"При анализе данных, полученных с использованием микрочипов, полученные измерения трактуют как непрерывные величины (лог-нормальное распределение). При анализе данных RNA-Seq, получаемые значения количества картируемых фрагментов натуральные, для анализа случайную величину принимают распределенной по Пуассону, как обратное биномиальное и даже бета-биномиальное."
