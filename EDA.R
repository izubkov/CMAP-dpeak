library(tidyverse)
library(magrittr)
library(plotly)
library(purrr)
library(tictoc)

# load all ----------------------------------------------------------------

DPK.files <-
  list.files("competitor_pack_v2/input/DPK.CP001_A549_24H_X1_B42/",
             full.names = T)
LIT.files <-
  list.files("competitor_pack_v2/input/LITMUS.KD017_A549_96H_X1_B42/",
             full.names = T)
DPK.all.txt <-
  map_dfr(DPK.files, read_tsv)
LIT.all.txt <-
  map_dfr(LIT.files, read_tsv)

all(DPK.all.txt %>%
      select(barcode_id) %>%
      distinct() ==
      LIT.all.txt %>%
      select(barcode_id) %>%
      distinct())

barcode_to_gene_map.txt <-
  read_tsv("submission/barcode_to_gene_map.txt", col_types = "iii")

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

source("scoring.R")

report(DPK.DECONV.uni, DPK.DECONV.out, DPK.DE.uni, DPK.DE.out)

# COR, AUC visualization --------------------------------------------------

cors <- compute_spearman_accuracy(DPK.DECONV.uni, DPK.DECONV.out)
cors <- compute_spearman_accuracy(DPK.DECONV.uni, DPK.DECONV.matlab)

# 976 = 16 * 61
plot_ly(x = rep(1:61, 16), y = rep(61:1, times = 1, each = 61),
        z = cors, type = "heatmap", text = ~names(cors))

# hist. of individual and all ---------------------------------------------

hist_FI <- function(filenames = DPK.files, gct = DPK.gct, out = DPK.out.gct,
                    txt = "A03") {
  cat("Barcodes of", txt, "\n")

  fid <- which(grepl(txt, filenames))
  stopifnot(length(fid) == 1)

  d <-
    read_tsv(filenames[fid],
             col_types = "dd",
             col_names = c("barcode_id", "FI"), skip = 1)

  li <-
    d %>%
    distinct(barcode_id) %>%
    select(barcode_id) %>%
    unlist(use.names = F)
  print(li)

  bid <- readline("Barcode: ")

  d %>%
    filter(barcode_id == bid) %>%
    select(FI) %>%
    unlist(use.names = F) %>%
    hist(breaks = length(.),
         main = str_c(bid, ", ", txt, ".txt"))

  genes <- barcode_to_gene_map.txt %>%
    filter(barcode_id == bid) %>% arrange(desc(high_prop))
  print(genes)

  gene_hi <- genes$gene_id[1]
  gene_lo <- genes$gene_id[2]

  ground_hi <- gct[gct$id == gene_hi, ][txt][[1]]
  ground_lo <- gct[gct$id == gene_lo, ][txt][[1]]
  out_hi <- out[out$id == gene_hi, ][txt][[1]]
  out_lo <- out[out$id == gene_lo, ][txt][[1]]
  abline(v = ground_hi, col = "red")
  abline(v = ground_lo, col = "blue")
  abline(v = out_hi, col = "darkred", lty = "dotted")
  abline(v = out_lo, col = "darkblue", lty = "dotted")
  cat("Ground hi/lo:", ground_hi, ground_lo, "\n")
  cat("Matlab out hi/lo:", out_hi, out_lo)
}

hist_FI_all <- function(d.all = DPK.all.txt, bid = 41) {
  d.all %>%
    filter(barcode_id == bid) %>%
    plot_ly(x = ~FI, type = "histogram")
}

# clusters ----------------------------------------------------------------

library(Gmedian)

kmeans_2_1_sizes <- function(d.all) {

  calc_kmeans <- function(x) {
    sizes <- kmeans(x, 2, algorithm = "Lloyd") %>% .[["size"]]
    if(sizes[1] > sizes[2]) {
      sizes[1] / sizes[2]
    } else {
      sizes[2] / sizes[1]
    }
  }

  tic("kmeans_2_1_sizes")

  sizes_2_1 <- d.all %>%
    group_by(barcode_id) %>%
    #filter(FI > 0.005 * max(FI)) %>%
    summarise(ratio = calc_kmeans(FI))

  toc()

  sizes_2_1
}

ratios <- kmeans_2_1_sizes(DPK.all.txt)
ratios %>%
  filter(ratio < 1.9 | ratio > 2.5) %T>%
  {.[["barcode_id"]] ->> outliers} %T>%
  {.[["ratio"]] %>% plot()} %>%
  dim()

for(i in seq_along(1:length(outliers))) {
  print(hist_FI_all(DPK.all.txt, outliers[i]))
}
