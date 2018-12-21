library(tidyverse)
library(purrr)

DPK.files <- list.files("competitor_pack_v2/input/DPK.CP001_A549_24H_X1_B42/",
                        full.names = T)
DPK.txt <-
  map_dfr(DPK.files, read_tsv)

LIT.files <- list.files("competitor_pack_v2/input/LITMUS.KD017_A549_96H_X1_B42/",
                        full.names = T)
LIT.txt <-
  map_dfr(DPK.files, read_tsv)

all(DPK.txt %>%
      select(barcode_id) %>%
      distinct() ==
      LIT.txt %>%
      select(barcode_id) %>%
      distinct())

barcode_to_gene_map.txt <-
  read_tsv("competitor_pack_v2/input/barcode_to_gene_map.txt")
