library(Seurat)



analysis_id <- paste0(
  "prepare_azimuth_reference_v1"
)
out_dir <- paste0("results/", analysis_id, "_results")
out_int_dir <- paste0("results/", analysis_id, "_intermediate_files")

if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

if (!dir.exists(out_int_dir)) {
  dir.create(out_int_dir)
}

pp_seu_list <- readRDS(
  "results/preprocess_v2_results/pp_seu_list.rds")

bbi <- pp_seu_list[["co"]]

# code adapted from
# <https://github.com/satijalab/azimuth-references/blob/b8b07dcdfcc09816a85aad07362e7bad4de03976/human_fetus/scripts/preprocess_BBI_reference.R>
Idents(object = bbi) <- "Organ_cell_lineage"

# normalize
bbi <- SCTransform(
  object = bbi,
  method = "glmGamPoi", 
  clip.range = c(-10000, 10),
  do.correct.umi = FALSE,
  conserve.memory = TRUE
)

# dimensional reduction
bbi <- RunPCA(object = bbi, npcs = 100)
bbi <- RunUMAP(object = bbi, dims = 1:100, return.model = TRUE)

saveRDS(bbi, file.path(out_int_dir, "intermediate_bbi.rds"))
