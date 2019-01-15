library(tidyverse)
library(magrittr)
library(plotly)
library(purrr)
library(tictoc)

source("R/functions.R")
source("R/scoring.R")

# load all ----------------------------------------------------------------

source("R/load.R")

# ground truth ------------------------------------------------------------

DPK.DE.uni     <- parse.gctx("ground-truth/DPK.CP001_A549_24H_X1_B42_DE_UNI.gct")
DPK.DECONV.uni <- parse.gctx("ground-truth/DPK.CP001_A549_24H_X1_B42_DECONV_UNI.gct")

LIT.DE.uni     <- parse.gctx("ground-truth/LITMUS.KD017_A549_96H_X1_B42_DE_UNI.gct")
LIT.DECONV.uni <- parse.gctx("ground-truth/LITMUS.KD017_A549_96H_X1_B42_DECONV_UNI.gct")

# matlab ------------------------------------------------------------------

DPK.DE.matlab     <- parse.gctx("competitor_pack_v2/output/DPK.CP001_A549_24H_X1_B42/level4_ZSPC_n376x976.gct")
DPK.DECONV.matlab <- parse.gctx("competitor_pack_v2/output/DPK.CP001_A549_24H_X1_B42.gct")

LIT.DE.matlab     <- parse.gctx("competitor_pack_v2/output/LITMUS.KD017_A549_96H_X1_B42/level4_ZSPC_n374x976.gct")
LIT.DECONV.matlab <- parse.gctx("competitor_pack_v2/output/LITMUS.KD017_A549_96H_X1_B42.gct")

# output ------------------------------------------------------------------

DPK.DE.out     <- parse.gctx("output/DPK.CP001_A549_24H_X1_B42/level4_ZSPC_n376x976.gct")
DPK.DECONV.out <- parse.gctx("output/DPK.CP001_A549_24H_X1_B42.gct")

LIT.DE.out     <- parse.gctx("output/LITMUS.KD017_A549_96H_X1_B42/level4_ZSPC_n374x976.gct")
LIT.DECONV.out <- parse.gctx("output/LITMUS.KD017_A549_96H_X1_B42.gct")

# COR, AUC reports --------------------------------------------------------

report(DPK.DECONV.uni, DPK.DECONV.out, DPK.DE.uni, DPK.DE.out)

# COR, AUC visualization --------------------------------------------------

COR <- compute_spearman_accuracy(DPK.DECONV.uni, DPK.DECONV.out)
COR <- compute_spearman_accuracy(DPK.DECONV.uni, DPK.DECONV.matlab)

# 976 = 16 * 61
plot_ly(x = rep(1:61, 16), y = rep(61:1, times = 1, each = 61),
        z = COR, type = "heatmap", text = ~names(COR))

# histograms --------------------------------------------------------------


# DEBUG stats -------------------------------------------------------------

load("dstat.RData")

# overview
dstat %>%
  group_by(reason) %>%
  count() %>%
  arrange(desc(n))

# persentage of bad clusters (max = 39.1 %)
dstat %>%
  group_by(bid) %>%
  summarize(
    bad_clusters =
      sum(reason %in% c("is.null(k)", "length(FI.n) < 2")) / n * 100) %>%
  arrange(desc(bad_clusters))

# read/write 100% accuracy ------------------------------------------------

#cmapR::parse.gctx(str_c("ground-truth/", args[1], "_DECONV_UNI.gct")) %>%
#  cmapR::write.gct(str_c(args[2], "/", args[1], ".gct"), appenddim = F)
