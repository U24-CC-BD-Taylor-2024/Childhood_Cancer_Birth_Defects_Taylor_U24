library(Seurat)
library(tidyverse)



analysis_id <- paste0(
  "azimuth_v2"
)
out_dir <- paste0("results/", analysis_id, "_results")
out_int_dir <- paste0("results/", analysis_id, "_intermediate_files")

if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

if (!dir.exists(out_int_dir)) {
  dir.create(out_int_dir)
}

options(future.globals.maxSize = 10 * 1024 ^ 3)

pp_seu_list <- readRDS(
  "results/preprocess_v2_results/pp_seu_list.rds")


map(pp_seu_list, dim)
map(pp_seu_list, function(xseu) {
  xcm <- GetAssayData(xseu[["RNA"]], layer = "counts")
  return(sum(rowSums(xcm > 0) > 0))
})
table(pp_seu_list$kd$stage)

run_az_seu_list <- pp_seu_list[
  c("kd", "cz", "ag_scRNA")]

walk(run_az_seu_list, function(x) {
  stopifnot(!is.null(x))
})

az_list <- imap(
  run_az_seu_list,
  function(xseu, xname) {
    print(xname)
    res <- Azimuth::RunAzimuth(
      xseu,
      file.path(
        "results/prepare_azimuth_reference_v1_results",
        "fetus_heart_ref"))

    cat("------------------\n\n\n\n\n\n")
    return(res)
  }
)


saveRDS(az_list, file.path(out_int_dir, "az_list.rds"))
