#!/usr/bin/env Rscript

# installed with cmapR
library(dplyr)
library(cmapR)

# installed with docker
library(doParallel)
library(foreach)
library(getopt)
library(magrittr)
library(purrr)
library(readr)
library(stringr)

source("R/kmeanspp.R")

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

DATA.all.txt <-
  map_dfr(DATA.files, read_tsv)

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

  FI <- FI[FI > 10] # filter lower threshold

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
    hi <- lo <- -median(FI)
    dstat[nrow(dstat)+1,] <<- c(bid, plate_name, "is.null(k)")
  } else {
    FI.1 <- FI[k$cluster == 1]
    FI.2 <- FI[k$cluster == 2]

    if(length(FI.1) < 2 || length(FI.2) < 2) {
      hi <- lo <- -median(FI)
      dstat[nrow(dstat)+1,] <<- c(bid, plate_name, "length(FI.n) < 2")
    } else {
      beads.1 <- length(FI.1)
      beads.2 <- length(FI.2)

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

save(dstat, file = "dstat.RData")
