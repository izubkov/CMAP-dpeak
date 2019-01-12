# barcode ids -------------------------------------------------------------

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

# hist. of individual and all ---------------------------------------------

hist_FI <- function(filenames = DPK.files,
                    ground.gct = DPK.DECONV.uni,
                    out = DPK.DECONV.out,
                    txt = "A03") {
  cat("Barcodes of", txt, "\n")

  fid <- which(grepl(txt, filenames))
  stopifnot(length(fid) == 1)

  d <-
    read_tsv(filenames[fid],
             col_types = "ii",
             col_names = c("barcode_id", "FI"), skip = 1)

  li <-
    d %>%
    distinct(barcode_id) %>%
    select(barcode_id) %>%
    unlist(use.names = F)
  print(li)

  bid <- readline("Barcode: ")

  # histogram
  d %>%
    filter(barcode_id == bid) %>%
    select(FI) %>%
    unlist(use.names = F) %>%
    hist(breaks = length(.),
         main = str_c(bid, ", ", txt, ".txt"))

  genes <- genes_by_bid(bid)
  gene_hi <- genes$gene_id[1] %>% as.character()
  gene_lo <- genes$gene_id[2] %>% as.character()
  print(genes)

  ground_hi <- ground.gct@mat[gene_hi, txt]
  ground_lo <- ground.gct@mat[gene_lo, txt]
  out_hi <- out@mat[gene_hi, txt]
  out_lo <- out@mat[gene_lo, txt]
  abline(v = ground_hi, col = "red")
  abline(v = ground_lo, col = "blue")
  abline(v = out_hi, col = "darkred", lty = "dotted")
  abline(v = out_lo, col = "darkblue", lty = "dotted")
  cat("Ground hi/lo:", ground_hi, ground_lo, "\n")
  cat("Matlab out hi/lo:", out_hi, out_lo)
}

hist_FI_all <- function(d.all = DPK.all.txt, bid = 12) {
  d.all %>%
    filter(barcode_id == bid) %>%
    plot_ly(x = ~FI, type = "histogram")
}
