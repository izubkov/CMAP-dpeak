#!/usr/bin/env Rscript

#
# TODO: foreach for parallel execution
# TODO: Gmeadian
# TODO: higher peak means higher proportion
#

# installed with cmapR
library(dplyr)
library(cmapR)

# installed with submission
library(doParallel)
library(foreach)
library(getopt)
library(magrittr)
library(purrr)
library(readr)
library(stringr)

# do not import into final submission
#library(tidyverse)
#library(tictoc)
#library(pryr)

args_spec <- matrix(c(
  "dspath",        "x", 1, "character",
  "out",           "y", 1, "character",
  "create_subdir", "z", 2, "integer",
  "plate",         "f", 2, "character"
), byrow = T, ncol = 4)
opt <- getopt::getopt(args_spec)

print(opt$dspath)
print(opt$out)

# load all ----------------------------------------------------------------

DATA.files <-
  list.files(str_c(opt$dspath, "/"),
             full.names = T)

extract_plate_name <- function(ss) {
  strsplit(ss, "_") %>%
    unlist() %>%
    tail(1) %>%
    strsplit(".txt") %>%
    unlist()
}

DATA.plates <-
  {
    DATA.files %>%
    lapply(extract_plate_name) %>%
    unlist()
  }

# DATA.all.txt <-
#   map_dfr(DATA.files, read_tsv)

barcode_to_gene_map.txt <-
  read_tsv("barcode_to_gene_map.txt", col_types = "iii")

# submission --------------------------------------------------------------

# TODO: just hardcode 11, 499 barcodes?
barcodes_to_skip <-
  barcode_to_gene_map.txt %>%
  group_by(barcode_id) %>%
  count() %T>%
  { filter(., n == 2) %>%
      .[["barcode_id"]] ->> barcodes } %>%
  filter(n < 2) %>%
  .[["barcode_id"]]

run_alg <- function(bid, FI) {
# fc <- new("flexclustControl", iter.max = 10, verbose = 1,
#           initcent = "kmeanspp")
# cc <- cclust(FI, 2, dist = "euclidean", method = "kmeans",
#              control = fc)

distEuclidean <- function(x, centers)
{
  if(ncol(x)!=ncol(centers))
    stop(sQuote("x")," and ",sQuote("centers"),
         " must have the same number of columns")
  z <- matrix(0, nrow=nrow(x), ncol=nrow(centers))
  for(k in 1:nrow(centers)){
    z[,k] <- sqrt( colSums((t(x) - centers[k,])^2) )
  }
  z
}

# TODO: sample from peaks/maximum values?
# TODO: different distance (manhattan)
kmeanspp <- function(x, k, family = NULL)
{
  centers <- matrix(0, nrow=k, ncol=ncol(x))
  centers[1,] <- x[sample(1:nrow(x), 1), , drop=FALSE]
  d <- distEuclidean(x, centers[1L,,drop=FALSE])^2
  for(l in 2:k){
    centers[l,] <- x[sample(1:nrow(x), 1, prob=d), , drop=FALSE]
    d <- pmin(d, distEuclidean(x, centers[l,,drop=FALSE])^2)
  }
  centers
}

  k <-
    tryCatch(
      {
        set.seed(42)
        cs <- kmeanspp(as.matrix(FI), 2)
        kmeans(x = FI, centers = cs, algorithm = "Lloyd")
      },
      warning = function(w) {
        NULL
      })

  if(is.null(k)) {
    hi <- lo <- median(FI)
  } else {
    x.1 <- FI[k$cluster == 1]
    x.2 <- FI[k$cluster == 1]
    max.1 <- max(x.1)
    max.2 <- max(x.2)

    if(max.1 > max.2) {
      hi <- median(FI[k$cluster == 1])
      lo <- median(FI[k$cluster == 2])
    } else {
      hi <- median(FI[k$cluster == 2])
      lo <- median(FI[k$cluster == 1])
    }
  }

  genes <-
    barcode_to_gene_map.txt %>%
    filter(barcode_id == bid)

  gene_hi <-
    genes[genes$high_prop == 1, ] %>%
    .[["gene_id"]] %>%
    as.character()

  gene_lo <-
    genes[genes$high_prop == 0, ] %>%
    .[["gene_id"]] %>%
    as.character()

  li <- vector(mode = "list", length = 2)
  names(li) <- c(gene_hi, gene_lo)
  li[[gene_hi]] <- hi
  li[[gene_lo]] <- lo

  return(li)
}

single_plate_processing <- function(filename) {
  plate_name <- extract_plate_name(filename)

  d <- read_tsv(filename, col_types = "ii")
  xs <-
    d %>%
    group_by(barcode_id) %>%
    filter(!barcode_id %in% barcodes_to_skip) %>%
    arrange(barcode_id) %>%
    split(.$barcode_id)

  col <-
    (foreach(x = xs) %do%
      run_alg(x$barcode_id[1], x$FI)) %>%
    unlist()

  return(col)
}

registerDoParallel(cores = detectCores(all.tests = T))
cols <-
  foreach(filename = DATA.files) %dopar%
    single_plate_processing(filename)
mat <- matrix(cols %>% unlist(),
              nrow = length(cols[[1]]),
              ncol = length(DATA.plates))
rownames(mat) <- names(cols[[1]])
colnames(mat) <- DATA.plates

# save DATA ---------------------------------------------------------------

output_name <-
  opt$dspath %>%
  strsplit("/") %>%
  unlist() %>%
  last()

gct <- new("GCT", mat=mat)
print("Saving GCT...")
gct %>%
  cmapR::write.gct(str_c(opt$out, "/", output_name, ".gct"), appenddim = F)

# all plates --------------------------------------------------------------

plates_all_processing <- function(d.all, debug = F) {

  run_kmeans <- function(f, li) {
    k <- do.call(what = f, args = li)
  }

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
#sol <- plates_all_processing(DATA.all.txt, debug = T)

#temp_gct <-
#  cmapR::parse.gctx(str_c("ground-truth/", args[1], "_DECONV_UNI.gct"))
#
#m <- ncol(temp_gct@mat)
#for(r in rownames(temp_gct@mat)) {
#  val <- sol[r][[1]]
#  temp_gct@mat[r,] <- rep(val, m) + rnorm(m)
#}
#
#print("Saving GCT...")
#temp_gct %>%
#  cmapR::write.gct(str_c(args[2], "/", args[1], ".gct"), appenddim = F)


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
