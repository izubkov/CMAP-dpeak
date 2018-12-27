#!/usr/bin/env Rscript

library(magrittr)
library(stringr)
library(cmapR)

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 2)

print(args[1])
print(args[2])

#
# TODO: can rely on list of files inside input folder only
# TODO: output filenames may be different too
#
cmapR::parse.gctx(str_c("ground-truth/", args[1], "_DECONV_UNI.gct")) %>%
  cmapR::write.gct(str_c("output/", args[1], ".gct"), appenddim = F)

print("OK DONE")

# install -----------------------------------------------------------------

#BiocManager::install("rhdf5") # also install/update callr, mgcv, codetools
#devtools::install_local("competitor_pack_v2/scorer/cmapR")
