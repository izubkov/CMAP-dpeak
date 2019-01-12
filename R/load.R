# load all ----------------------------------------------------------------

DPK.files <-
  list.files("competitor_pack_v2/input/DPK.CP001_A549_24H_X1_B42/",
             full.names = T)

DPK.all.txt <-
  map_dfr(DPK.files, read_tsv, col_types = "ii")

barcode_to_gene_map.txt <-
  read_tsv("submission/barcode_to_gene_map.txt", col_types = "iii")

# LIT.all.txt <-
#   map_dfr(LIT.files, read_tsv)
#
# LIT.files <-
#   list.files("competitor_pack_v2/input/LITMUS.KD017_A549_96H_X1_B42/",
#              full.names = T)
#
# all(DPK.all.txt %>%
#       select(barcode_id) %>%
#       distinct() ==
#       LIT.all.txt %>%
#       select(barcode_id) %>%
#       distinct())
