#!/usr/bin/env Rscript

#
# TODO: can rely on list of files inside input folder only
# TODO: output filenames may be different too
#
# TODO: foreach for parallel execution
#

library(magrittr)
library(stringr)
library(cmapR)

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 2)

print(args[1])
print(args[2])

plate_wide_processing <- function(d.all) {
  # TODO: matrix?
  res <- tibble(gene_hi = character(), center_hi = double(),
                gene_lo = character(), center_lo = double())

  bs <-
    d.all %>%
    distinct(barcode_id) %>%
    .[["barcode_id"]]

  for(i in seq_along(1:length(bs))) {
    barcode <- bs[i]
    xs <-
      d.all %>%
      filter(barcode_id == barcode) %>%
      .[["FI"]]

    li <- list(x = xs, centers = 2, algorithm = "MacQueen")
    k <- run_kmeans(kmeans, li)

    # ...
  }
}

run_kmeans <- function(f, li) {
  k <- do.call(what = f, args = li)
}

k <- kmeans(xs, 2, algorithm = "MacQueen")
k$centers # [1] [2]
k$size # [1] [2]

barcode_to_gene_map.txt %>% filter(barcode_id == 64)

temp_gct <-
  cmapR::parse.gctx("ground-truth/DPK.CP001_A549_24H_X1_B42_DECONV_UNI.gct")

cmapR::parse.gctx(str_c("ground-truth/", args[1], "_DECONV_UNI.gct")) %>%
  cmapR::write.gct(str_c("output/", args[1], ".gct"), appenddim = F)

# install -----------------------------------------------------------------

#BiocManager::install("rhdf5") # also install/update callr, mgcv, codetools
#devtools::install_local("competitor_pack_v2/scorer/cmapR")
