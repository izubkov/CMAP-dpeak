#! /usr/bin/env Rscript

library(magrittr)
library(cmapR)

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 2)

print(args[1])
print(args[2])

cmapR::parse.gctx("ground-truth/DPK.CP001_A549_24H_X1_B42_DECONV_UNI.gct") %>%
  cmapR::write.gct("output/DPK.CP001_A549_24H_X1_B42.gct", appenddim = F)
cmapR::parse.gctx("ground-truth/LITMUS.KD017_A549_96H_X1_B42_DECONV_UNI.gct") %>%
  cmapR::write.gct("output/LITMUS.KD017_A549_96H_X1_B42.gct", appenddim = F)

# install -----------------------------------------------------------------

#BiocManager::install("rhdf5") # also install/update callr, mgcv, codetools
#devtools::install_local("competitor_pack_v2/scorer/cmapR")
