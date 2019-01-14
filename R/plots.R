require(plotly)

# plot density ------------------------------------------------------------

plot_densities <- function(FI,
                           fun = function(x) x,
                           kmeans.result = k) {
  FI.t <- fun(FI)
  FI.1 <- FI.t[kmeans.result$cluster == 1]
  FI.2 <- FI.t[kmeans.result$cluster == 2]
  ds.0 <- density(FI.t, bw = "SJ", kernel = "biweight", adjust = 1)
  ds.1 <- density(FI.1, bw = "SJ", kernel = "biweight", adjust = 1)
  ds.2 <- density(FI.2, bw = "SJ", kernel = "biweight", adjust = 1)
  max_y <- max(c(ds.0$y, ds.1$y, ds.2$y))
  plot(ds.0, type = "n", ylim = c(0, max_y*1.1))
  lines(ds.0, col = "black")
  lines(ds.1, col = "green")
  lines(ds.2, col = "red")
}
# plot_densities(FI)
# plot_densities(FI, log)

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
