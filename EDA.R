library(tidyverse)
library(purrr)

# load all ----------------------------------------------------------------

DPK.files <-
  list.files("competitor_pack_v2/input/DPK.CP001_A549_24H_X1_B42/",
             full.names = T)
LIT.files <-
  list.files("competitor_pack_v2/input/LITMUS.KD017_A549_96H_X1_B42/",
             full.names = T)
DPK.all.txt <-
  map_dfr(DPK.files, read_tsv)
LIT.all.txt <-
  map_dfr(DPK.files, read_tsv)

all(DPK.all.txt %>%
      select(barcode_id) %>%
      distinct() ==
      LIT.all.txt %>%
      select(barcode_id) %>%
      distinct())

barcode_to_gene_map.txt <-
  read_tsv("competitor_pack_v2/input/barcode_to_gene_map.txt")

# ground truth ------------------------------------------------------------

DPK.gct <- read_tsv("competitor_pack_v2/ground-truth/DPK.CP001_A549_24H_X1_B42_DECONV_UNI.gct", skip = 2)
DPK.DE.gct <- read_tsv("competitor_pack_v2/ground-truth/DPK.CP001_A549_24H_X1_B42_DE_UNI.gct", skip = 2)
LIT.gct <- read_tsv("competitor_pack_v2/ground-truth/LITMUS.KD017_A549_96H_X1_B42_DECONV_UNI.gct", skip = 2)
LIT.DE.gct <- read_tsv("competitor_pack_v2/ground-truth/LITMUS.KD017_A549_96H_X1_B42_DE_UNI.gct", skip = 2)

# output ------------------------------------------------------------------

DPK.out.gct <- read_tsv("competitor_pack_v2/output/DPK.CP001_A549_24H_X1_B42.gct", skip = 2)
LIT.out.gct <- read_tsv("competitor_pack_v2/output/LITMUS.KD017_A549_96H_X1_B42.gct", skip = 2)

read_FI <- function(file, list_barcodes = T) {
  cat("File", (file %>% str_split("_"))[[1]][c(7, 12)], "\n" )
  d <-
    read_tsv(file,
             col_types = "dd",
             col_names = c("barcode_id", "FI"), skip = 1)
  li <-
    d %>%
    distinct(barcode_id) %>%
    select(barcode_id) %>%
    unlist(use.names = F)
  if(list_barcodes) print(li)
  return(d)
}

hist_FI <- function(d, bid) {
  ids <- d %>% distinct(barcode_id) %>% unlist(use.names = F)
  d %>%
    filter(barcode_id == ids[bid]) %>%
    select(FI) %>%
    unlist(use.names = F) %>%
    hist(breaks = length(.))
}

extract_name <- function(filename) {
  txt <- str_split(filename, "_")[[1]][12]
   %>%
    str_split(".")[[1]][1]
}

DPK.files[1] %>% read_FI() %>% hist_FI(490)

hist_FI <- function(filenames, gct, out, txt) {
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

  ground_hi <- gct[gct$id == gene_hi, ][txt][[1]]
  ground_lo <- gct[gct$id == gene_lo, ][txt][[1]]
  out_hi <- out[out$id == gene_hi, ][txt][[1]]
  out_lo <- out[out$id == gene_lo, ][txt][[1]]
  abline(v = ground_hi, col = "red")
  abline(v = ground_lo, col = "blue")
  abline(v = out_hi, col = "darkred")
  abline(v = out_lo, col = "darkblue")
  cat("Ground hi/lo:", ground_hi, ground_lo, "\n")
  cat("Matlab out hi/lo:", out_hi, out_lo)
}
