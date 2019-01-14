# barcode ids -------------------------------------------------------------

require(tidyverse)

bids_by_gene <- function(genes, d = barcode_to_gene_map.txt) {
  d[d$gene_id %in% as.integer(genes), ] %>%
    arrange(desc(high_prop))
}

genes_by_bid <- function(bids, d = barcode_to_gene_map.txt) {
  d[d$barcode_id %in% as.integer(bids), ] %>%
    arrange(desc(high_prop))
}

intensities <- function(filenames = DPK.files, bid = 12, txt = "A11") {
  fid <- which(grepl(txt, filenames))
  stopifnot(length(fid) == 1)

  d <-
    read_tsv(filenames[fid],
             col_types = "ii",
             col_names = c("barcode_id", "FI"), skip = 1) %>%
    filter(barcode_id == bid)

  d$FI
}

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

# install -----------------------------------------------------------------

#BiocManager::install("rhdf5") # also install/update callr, mgcv, codetools
#devtools::install_local("competitor_pack_v2/scorer/cmapR")

# clear warnings
#assign("last.warning", NULL, envir = baseenv())

# stop on warnings
#options(warn = 2)
