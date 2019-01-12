# barcode ids -------------------------------------------------------------

bids_by_gene <- function(genes, d = barcode_to_gene_map.txt) {
  d[d$gene_id %in% as.integer(genes), ]
}

genes_by_bid <- function(bids, d = barcode_to_gene_map.txt) {
  d[d$barcode_id %in% as.integer(bids), ]
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

  ground_hi <- ground.gct[ground.gct$id == gene_hi, ][txt][[1]]
  ground_lo <- ground.gct[ground.gct$id == gene_lo, ][txt][[1]]
  out_hi <- out[out$id == gene_hi, ][txt][[1]]
  out_lo <- out[out$id == gene_lo, ][txt][[1]]
  abline(v = ground_hi, col = "red")
  abline(v = ground_lo, col = "blue")
  abline(v = out_hi, col = "darkred", lty = "dotted")
  abline(v = out_lo, col = "darkblue", lty = "dotted")
  cat("Ground hi/lo:", ground_hi, ground_lo, "\n")
  cat("Matlab out hi/lo:", out_hi, out_lo)
}

hist_FI_all <- function(d.all = DPK.all.txt, bid = 41) {
  d.all %>%
    filter(barcode_id == bid) %>%
    plot_ly(x = ~FI, type = "histogram")
}
