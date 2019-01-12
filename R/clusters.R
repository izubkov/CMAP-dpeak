# clusters ----------------------------------------------------------------

library(Gmedian) # TODO

kmeans_2_1_sizes <- function(d.all) {

  calc_kmeans <- function(x) {
    sizes <- kmeans(x, 2, algorithm = "Lloyd") %>% .[["size"]]
    if(sizes[1] > sizes[2]) {
      sizes[1] / sizes[2]
    } else {
      sizes[2] / sizes[1]
    }
  }

  tic("kmeans_2_1_sizes")

  sizes_2_1 <- d.all %>%
    group_by(barcode_id) %>%
    #filter(FI > 0.005 * max(FI)) %>%
    summarise(ratio = calc_kmeans(FI))

  toc()

  sizes_2_1
}

ratios <- kmeans_2_1_sizes(DPK.all.txt)
ratios %>%
  filter(ratio < 1.9 | ratio > 2.5) %T>%
  {.[["barcode_id"]] ->> outliers} %T>%
  {.[["ratio"]] %>% plot()} %>%
  dim()

for(i in seq_along(1:length(outliers))) {
  print(hist_FI_all(DPK.all.txt, outliers[i]))
}
