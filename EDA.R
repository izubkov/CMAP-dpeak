library(tidyverse)
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
  read_tsv("competitor_pack_v2/input/barcode_to_gene_map.txt")

# ground truth ------------------------------------------------------------

DPK.gct <- read_tsv("competitor_pack_v2/ground-truth/DPK.CP001_A549_24H_X1_B42_DECONV_UNI.gct", skip = 2)
DPK.DE.gct <- read_tsv("competitor_pack_v2/ground-truth/DPK.CP001_A549_24H_X1_B42_DE_UNI.gct", skip = 2)
LIT.gct <- read_tsv("competitor_pack_v2/ground-truth/LITMUS.KD017_A549_96H_X1_B42_DECONV_UNI.gct", skip = 2)
LIT.DE.gct <- read_tsv("competitor_pack_v2/ground-truth/LITMUS.KD017_A549_96H_X1_B42_DE_UNI.gct", skip = 2)

# output ------------------------------------------------------------------

DPK.out.gct <- read_tsv("competitor_pack_v2/output/DPK.CP001_A549_24H_X1_B42.gct", skip = 2)
LIT.out.gct <- read_tsv("competitor_pack_v2/output/LITMUS.KD017_A549_96H_X1_B42.gct", skip = 2)

# hist. of individual and all ---------------------------------------------

hist_FI <- function(filenames, gct, out, txt) {
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

hist_FI_all <- function(d.all, bid) {
  d.all %>%
    filter(barcode_id == bid) %>%
    plot_ly(x = ~FI, type = "histogram")
}

# clusters ----------------------------------------------------------------

library(Gmedian)

kmeans_2_1_sizes <- function(d.all) {
  tic("kmeans_2_1_sizes")
  calc_kmeans <- function(x) {
    sizes <- kmeans(x, 2, algorithm = "MacQueen") %>% .[["size"]]
    if(sizes[1] > sizes[2]) {
      sizes[1] / sizes[2]
    } else {
      sizes[2] / sizes[1]
    }
  }

  sizes_2_1 <- d.all %>%
    group_by(barcode_id) %>%
    summarise(calc_kmeans(FI)) -> result
  toc()
  sizes_2_1
}

# TODO: plot outliers
