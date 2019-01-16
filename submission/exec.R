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

  # TODO: filter upper threshold (?)
  FI <- FI[FI > 10] # filter lower threshold

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
  kmeanspp <- function(x, k = 2, start = "rand")
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
        cs <- kmeanspp(as.matrix(FI), 2)
        kmeans(x = FI, centers = cs, algorithm = "Lloyd", iter.max = 15)
      },
      warning = function(w) {
        NULL
      })

  if(is.null(k)) {
    hi <- lo <- median(FI)
    #dstat[nrow(dstat)+1,] <<- c(bid, plate_name, "is.null(k)")
  } else {
    FI.1 <- FI[k$cluster == 1]
    FI.2 <- FI[k$cluster == 2]

    if(length(FI.1) < 2 || length(FI.2) < 2) {
      hi <- lo <- median(FI)
      #dstat[nrow(dstat)+1,] <<- c(bid, plate_name, "length(FI.n) < 2")
    } else {
      h.1 <- hist(FI.1, breaks = length(FI.1), plot = F)
      h.2 <- hist(FI.2, breaks = length(FI.2), plot = F)
      peak.1 <- h.1$breaks[h.1$counts == max(h.1$counts)][1]
      peak.2 <- h.2$breaks[h.2$counts == max(h.2$counts)][1]
      beads.1 <- length(FI.1)
      beads.2 <- length(FI.2)

      if(peak.1 > peak.2 && beads.1 > beads.2) {
        # sure the first cluster is high_prop
        hi <- median(FI[k$cluster == 1])
        lo <- median(FI[k$cluster == 2])
        #dstat[nrow(dstat)+1,] <<- c(bid, plate_name, "sure:1")
      } else if(peak.1 < peak.2 && beads.1 < beads.2) {
        # sure the second cluster is high_prop
        hi <- median(FI[k$cluster == 2])
        lo <- median(FI[k$cluster == 1])
        #dstat[nrow(dstat)+1,] <<- c(bid, plate_name, "sure:2")
      } else {
        # unsure case
        if(beads.1 > beads.2) {
          hi <- median(FI[k$cluster == 1])
          lo <- median(FI[k$cluster == 2])
          #dstat[nrow(dstat)+1,] <<- c(bid, plate_name, "unsure:1")
        } else {
          hi <- median(FI[k$cluster == 2])
          lo <- median(FI[k$cluster == 1])
          #dstat[nrow(dstat)+1,] <<- c(bid, plate_name, "unsure:2")
        }
      }

      # END OF: good clusters
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
    arrange(barcode_id) %>% # TODO: not arrange but return single column with names
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

# voting ------------------------------------------------------------------

genes_by_bids <- function(bids, d = barcode_to_gene_map.txt) {
  result <-
    d[d$barcode_id %in% as.integer(bids), ] %>%
    arrange(desc(high_prop)) %>%
    unlist()
  result[3:4] %>% as.character()
}

for(bid in barcodes) {
  genes <- genes_by_bids(bid, barcode_to_gene_map.txt)
  gene_hi <- genes[1]
  gene_lo <- genes[2]
  ms.hi <- mat[gene_hi,]
  ms.lo <- mat[gene_lo,]
  ms.hi <- ms.hi[ms.hi > 0]
  ms.lo <- ms.lo[ms.lo > 0]
  mean.hi <- median(ms.hi)
  mean.lo <- median(ms.lo)
  sd.hi <- sd(ms.hi)
  sd.lo <- sd(ms.lo)
  for(well in names(ms.hi)) {
    hi <- ms.hi[well]
    lo <- ms.lo[well]
    if(mean.lo - sd.lo <= hi && hi <= mean.lo + sd.lo &&
       mean.hi - sd.hi <= lo && lo <= mean.hi + sd.hi) {
      mat[gene_hi,well] <- lo
      mat[gene_lo,well] <- hi
    }
  }
}


# mat <- apply(t(mat), c(2, 1), function(x) if(x < 0) { -x } else x)

# for(bid in barcodes) {
#   genes <- genes_by_bids(bid, barcode_to_gene_map.txt)
#   gene_hi <- genes[1]
#   gene_lo <- genes[2]
#   ms.hi <- mat[gene_hi,]
#   ms.lo <- mat[gene_lo,]
#   ms.hi <- ms.hi[ms.hi > 0]
#   ms.lo <- ms.lo[ms.lo > 0]
#   mean.hi <- mean(ms.hi)
#   mean.lo <- mean(ms.lo)
#   for(well in names(ms.hi)) {
#     hi <- ms.hi[well]
#     if(hi < 0) {
#       mat[gene_hi,well] <- mean.hi
#       mat[gene_lo,well] <- mean.lo
#     }
#   }
# }

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
