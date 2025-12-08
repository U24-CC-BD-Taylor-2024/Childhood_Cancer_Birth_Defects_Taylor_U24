# This script is adapted from
# <https://github.com/satijalab/azimuth-references/blob/b8b07dcdfcc09816a85aad07362e7bad4de03976/human_fetus/scripts/make_BBI_reference.R>
library(Seurat)
library(stringr)
library(magrittr)
library(Azimuth)



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



print("Reading object")
bbi <- readRDS(file.path(out_int_dir, "intermediate_bbi.rds"))
print("Done reading in")


ocl_match <- str_match(
  bbi[['Organ_cell_lineage', drop=TRUE]],
  "^([^-]+)-(.+)$")

stopifnot(identical(0L, sum(is.na(ocl_match))))
annotation.l1 <- factor(ocl_match[, 3])

annotation.l2 <- factor(
  bbi[['Organ_cell_lineage', drop = TRUE]])

stopifnot(identical(
  unname(bbi[['Organ_cell_lineage', drop = TRUE]]),
  as.character(annotation.l2)))

stopifnot(identical(
  paste(
    bbi[['Organ', drop = TRUE]],
    as.character(annotation.l1),
    sep = "-"),
  as.character(annotation.l2)
))


organ <- factor(bbi[['Organ', drop = TRUE]])

stopifnot(identical(
  paste(
    as.character(organ),
    as.character(annotation.l1),
    sep = "-"),
  as.character(
    bbi[['Organ_cell_lineage', drop = TRUE]])
))

stopifnot(identical(0L, sum(is.na(organ))))
stopifnot(identical(0L, sum(is.na(annotation.l1))))
stopifnot(identical(0L, sum(is.na(annotation.l2))))

bbi[['annotation.l1']] <- annotation.l1
bbi[['annotation.l2']] <- annotation.l2

# subset plotref
set.seed(seed = 42)
plotting.cells <- sample(x = Cells(x = bbi), size = 10**5)
plotref <- subset(x = bbi[["umap"]], cells = plotting.cells)

# make azimuth reference
ref <- AzimuthReference(
    object = bbi,
    refUMAP = "umap",
    refDR = "pca",
    refAssay = "SCT",
    plotref = plotref,
    metadata = c('annotation.l1','annotation.l2'), 
    dims = 1:100,
    reference.version = "1.0.0"
)

# save
ref_out_dir <- file.path(out_dir, "fetus_heart_ref")
if (!dir.exists(ref_out_dir)) {
  dir.create(ref_out_dir)
}

saveRDS(
  object = ref,
  file = file.path(ref_out_dir, "ref.Rds"))

SaveAnnoyIndex(
  object = ref[['refdr.annoy.neighbors']],
  file = file.path(ref_out_dir, "idx.annoy"))

saveRDS(
  object = bbi,
  file = file.path(
    out_int_dir,
    "co_heart_azimuth_reference_seurat_object.rds"))
