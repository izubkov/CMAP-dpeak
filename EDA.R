library(tidyverse)
library(purrr)

DPK.files <- list.files("competitor_pack_v2/input/DPK.CP001_A549_24H_X1_B42/",
                        full.names = T)
DPK.all.txt <-
  map_dfr(DPK.files, read_tsv)

LIT.files <- list.files("competitor_pack_v2/input/LITMUS.KD017_A549_96H_X1_B42/",
                        full.names = T)
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

DPK.gct <-
  read_tsv("competitor_pack_v2/ground-truth/DPK.CP001_A549_24H_X1_B42_DECONV_UNI.gct", skip = 2)
