#!/usr/bin/env Rscript

#
# TODO: Gmeadian
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

# debug statistics
dstat <-
  data.frame(bid = integer(),
             plate_name = character(),
             reason = character(),
             stringsAsFactors = F)

run_alg <- function(bid, FI, plate_name = NULL) {
  # fc <- new("flexclustControl", iter.max = 10, verbose = 1,
  #           initcent = "kmeanspp")
  # cc <- cclust(FI, 2, dist = "euclidean", method = "kmeans",
  #              control = fc)

  FI <- FI[FI > 10] # filter low threshold

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

  # TODO: different distance (manhattan)?
  # TODO: try oversampling
  kmeanspp <- function(x, k = 2, start = "max")
  {
    centers <- matrix(0, nrow=k, ncol=ncol(x))
    if(start == "max") {
      h <- hist(x, breaks = length(x), plot = F)
      mx <- h$breaks[h$counts == max(h$counts)]
      centers[1,] <- sample(mx, 1)
    } else if(start == "rand") {
      centers[1,] <- x[sample(1:nrow(x), 1), , drop=FALSE]
    }
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
        cs <- kmeanspp(as.matrix(FI), 2, start = "max")
        kmeans(x = FI, centers = cs, algorithm = "Lloyd")
      },
      warning = function(w) {
        NULL
      })

  if(is.null(k)) {
    hi <- lo <- median(FI)
    dstat[nrow(dstat)+1,] <<- c(bid, plate_name, "is.null(k)")
  } else {
    FI.1 <- FI[k$cluster == 1]
    FI.2 <- FI[k$cluster == 2]

    if(length(FI.1) < 2 || length(FI.2) < 2) {
      hi <- lo <- median(FI)
      dstat[nrow(dstat)+1,] <<- c(bid, plate_name, "length(FI.n) < 2")
    } else {
      # TODO: remove (?)
      ds.1 <- ds.2 <- data.frame(y = 0)
      ds.1 <- try(density(FI.1, bw = "SJ", kernel = "gaussian", n = 64))
      ds.2 <- try(density(FI.2, bw = "SJ", kernel = "gaussian", n = 64))
      beads.1 <- length(FI.1)
      beads.2 <- length(FI.2)

      if(!is.atomic(ds.1) && !is.atomic(ds.2)) {
        peak.1 <- max(ds.1$y)
        peak.2 <- max(ds.2$y)

        if(peak.1 > peak.2 && beads.1 > beads.2) {
          # sure the first cluster is high_prop
          hi <- median(FI[k$cluster == 1])
          lo <- median(FI[k$cluster == 2])
          dstat[nrow(dstat)+1,] <<- c(bid, plate_name, "sure:1")
        } else if(peak.1 < peak.2 && beads.1 < beads.2) {
          # sure the second cluster is high_prop
          hi <- median(FI[k$cluster == 2])
          lo <- median(FI[k$cluster == 1])
          dstat[nrow(dstat)+1,] <<- c(bid, plate_name, "sure:2")
        } else {
          # unsure case
          if(beads.1 > beads.2) {
            hi <- median(FI[k$cluster == 1])
            lo <- median(FI[k$cluster == 2])
            dstat[nrow(dstat)+1,] <<- c(bid, plate_name, "unsure:1")
          } else {
            hi <- median(FI[k$cluster == 2])
            lo <- median(FI[k$cluster == 1])
            dstat[nrow(dstat)+1,] <<- c(bid, plate_name, "unsure:2")
          }
        }
      } else {
        # unsure case
        if(beads.1 > beads.2) {
          hi <- median(FI[k$cluster == 1])
          lo <- median(FI[k$cluster == 2])
          dstat[nrow(dstat)+1,] <<- c(bid, plate_name, "non_atomic:1")
        } else {
          hi <- median(FI[k$cluster == 2])
          lo <- median(FI[k$cluster == 1])
          dstat[nrow(dstat)+1,] <<- c(bid, plate_name, "non_atomic:2")
        }
      }
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
      run_alg(x$barcode_id[1], x$FI, plate_name)) %>%
    unlist()

  return(col)
}

registerDoParallel(cores = detectCores(all.tests = T))
cols <-
  # foreach(filename = DATA.files) %do%
  #   single_plate_processing(filename)
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

# save DEBUG --------------------------------------------------------------

#save(dstat, file = "dstat.RData")

# TODO --------------------------------------------------------------------

# # detect peaks
# FI.log2 <- log(FI, 2)
# # "epanechnikov", "rectangular", "triangular", "biweight", "cosine", "optcosine"
# ds <- density(FI.log2, bw = "SJ", adjust = 1,
#               kernel = "gaussian",
#               weights = NULL, window = kernel)
# plot(ds)
