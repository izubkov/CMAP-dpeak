#!/usr/bin/env Rscript

#
# TODO: foreach for parallel execution
#

library(magrittr)
library(purrr)
library(stringr)
library(cmapR)

# do not import into final solution
library(tidyverse)
library(tictoc)
library(pryr)

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 2)

print(args[1]) # input
print(args[2]) # output

# load all ----------------------------------------------------------------

DATA.files <-
  list.files(str_c("input/", args[1], "/"),
             full.names = T)

DATA.plates <-
  {
    DATA.files %>%
    lapply(function(x)
      strsplit(x, "_") %>%
      unlist() %>%
      tail(1) %>%
      strsplit(".txt") %>%
      unlist()) %>%
    unlist()
  }

DATA.all.txt <-
  map_dfr(DATA.files, read_tsv)

barcode_to_gene_map.txt <-
  read_tsv("input/barcode_to_gene_map.txt", col_types = "iii")

# solution ----------------------------------------------------------------

# TODO: just hardcode 11, 499 barcodes?
barcodes_to_skip <-
  barcode_to_gene_map.txt %>%
  group_by(barcode_id) %>%
  count() %T>%
  { filter(., n == 2) %>%
      .[["barcode_id"]] ->> barcodes } %>%
  filter(n < 2) %>%
  .[["barcode_id"]]

rows <-
  barcodes %>%
  lapply(function(x) {
    barcode_to_gene_map.txt[barcode_to_gene_map.txt$barcode_id == x, ] %>%
      .[["gene_id"]] %>%
      unlist()
    }) %>%
  unlist()
mat <- matrix(nrow = length(rows), ncol = length(DATA.plates))
colnames(mat) <- DATA.plates
rownames(mat) <- rows

single_plate_processing <- function(filename) {
  d <- read_tsv(filename, col_types = "ii")
  bids.all <-
    d %>%
    group_by(barcode_id) %>%
    nest()
  for(bid in bids.all$barcode_id) {
    xs <-
      bids.all[bids.all$barcode_id == bid, ]$data %>%
      unlist(use.names = F)
    print(xs)

    # TODO: ...
  }
}

run_kmeans <- function(f, li) {
  k <- do.call(what = f, args = li)
}

# TODO: medians of k-mean separated clusters
# TODO: individual medians of k-mean separated clusters
#
# TODO: Gmeadian
plates_all_processing <- function(d.all, debug = F) {

  res <- list()

  bs <-
    d.all %>%
    distinct(barcode_id) %>%
    .[["barcode_id"]]

  for(i in seq_along(1:length(bs))) {

    bid <- bs[i]

    xs <-
      d.all %>%
      filter(barcode_id == bid) %>%
      .[["FI"]]

    li <- list(x = xs, centers = 2, algorithm = "MacQueen")

    k <-
      tryCatch(
        { run_kmeans(kmeans, li) },
        warning = function(w) {
          # When kmeans does not converges it prints warning and
          # does not returns anything so k still points to prev. value.
          # Use plate-wide median.
          if(debug) {
            cat("Closure", address(k), "\n")
            cat(bid, "does not converges in 10 iterations\n")
          }
          med <- median(xs)
          list(size = -1, centers = med)
        })

    # assign clusters
    if(k$size == -1) {
      hi <- lo <- k$centers
    } else {
      x.1 <- xs[k$cluster == 1]
      x.2 <- xs[k$cluster == 1]
      max.1 <- max(x.1)
      max.2 <- max(x.2)

      # TODO: higher peak means higher proportion
      if(max.1 > max.2) {
        hi <- median(xs[k$cluster == 1])
        lo <- median(xs[k$cluster == 2])
      } else {
        hi <- median(xs[k$cluster == 2])
        lo <- median(xs[k$cluster == 1])
      }
    }

    genes <- barcode_to_gene_map.txt %>% filter(barcode_id == bid)
    if(nrow(genes) == 2) { # exclude barcodes 11 and 499
      gene.hi <-
        genes[genes$high_prop == 1, ] %>%
        .[["gene_id"]] %>%
        as.character()
      gene.lo <-
        genes[genes$high_prop == 0, ] %>%
        .[["gene_id"]] %>%
        as.character()

      res[gene.hi] <- hi
      res[gene.lo] <- lo
    }
  }
  res
}

sol <- plates_all_processing(DATA.all.txt, debug = T)

# save DATA ---------------------------------------------------------------

mat <- matrix(stats::rnorm(100), ncol=10)
rownames(mat) <- letters[1:10]
colnames(mat) <- LETTERS[1:10]
(my_ds <- new("GCT", mat=mat))

temp_gct <-
  cmapR::parse.gctx(str_c("ground-truth/", args[1], "_DECONV_UNI.gct"))

m <- ncol(temp_gct@mat)
for(r in rownames(temp_gct@mat)) {
  val <- sol[r][[1]]
  temp_gct@mat[r,] <- rep(val, m) + rnorm(m)
}

print("Saving GCT...")
temp_gct %>%
  cmapR::write.gct(str_c(args[2], "/", args[1], ".gct"), appenddim = F)

# read/write 100% accuracy ------------------------------------------------

#cmapR::parse.gctx(str_c("ground-truth/", args[1], "_DECONV_UNI.gct")) %>%
#  cmapR::write.gct(str_c(args[2], "/", args[1], ".gct"), appenddim = F)

# install -----------------------------------------------------------------

#BiocManager::install("rhdf5") # also install/update callr, mgcv, codetools
#devtools::install_local("competitor_pack_v2/scorer/cmapR")

# clear warnings
#assign("last.warning", NULL, envir = baseenv())

# stop on warnings
#options(warn = 2)
