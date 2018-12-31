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

# TODO: medians of k-mean separated clusters
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
    #k <- tryCatch({run_kmeans(kmeans, li)},
    #              warning = function(w) {
    #                cat(barcode, "does not converges in 10 iterations")
    #                data.frame(size = c(0, 0), centers = c(0, 0))
    #              })
    if(k$size[1] > k$size[2]) {
      hi <- k$centers[1]
      lo <- k$centers[2]
    } else {
      hi <- k$centers[2]
      lo <- k$centers[1]
    }

    genes <- barcode_to_gene_map.txt %>% filter(barcode_id == barcode)
    if(nrow(genes) == 2) { # barcodes 11 and 499
      gene.hi <-
        genes[genes$high_prop == 1, ] %>%
        .[["gene_id"]] %>%
        as.character()
      gene.lo <-
        genes[genes$high_prop == 0, ] %>%
        .[["gene_id"]] %>%
        as.character()

      res %<>% add_row(gene_hi = gene.hi, center_hi = hi,
                       gene_lo = gene.lo, center_lo = lo)
    }
  }
  res
}

run_kmeans <- function(f, li) {
  k <- do.call(what = f, args = li)
}

tic("DPK")
dpk <- plate_wide_processing(DPK.all.txt)
toc()
tic("LIT")
lit <- plate_wide_processing(LIT.all.txt)
toc()

temp_gct <-
  cmapR::parse.gctx("ground-truth/DPK.CP001_A549_24H_X1_B42_DECONV_UNI.gct")

cmapR::parse.gctx(str_c("ground-truth/", args[1], "_DECONV_UNI.gct")) %>%
  cmapR::write.gct(str_c("output/", args[1], ".gct"), appenddim = F)

# install -----------------------------------------------------------------

#BiocManager::install("rhdf5") # also install/update callr, mgcv, codetools
#devtools::install_local("competitor_pack_v2/scorer/cmapR")

# clear warnings
assign("last.warning", NULL, envir = baseenv())

# stop on warnings
options(warn = 2)
