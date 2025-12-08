library(Seurat)
library(tidyverse)



analysis_id <- paste0(
  "preprocess_v2"
)
out_dir <- paste0("results/", analysis_id, "_results")

if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}



kd_rds <- readRDS(
  file.path(
    "data/scrna_heart_development/GSE216019_Knight_Davaapil",
    "GSE216019_Adult_Foetal_Integrated_Hearts_VKS2022.rds"))

cz_tbl <- read_tsv(
  file.path(
    "data/scrna_heart_development/GSE106118_Cui_Zheng",
    "GSE106118_UMI_count_merge.txt"))

ag_data_dir <- file.path(
  "data/scrna_heart_development",
  "EGAS00001003996_asp_giacomello",
  "dev_heart_count_matrices_and_meta/Filtered"
)

ag_list <- list(
  scRNA = list(
    count_matrix = read_tsv(
      file.path(
        ag_data_dir,
        "Developmental_heart_filtered_scRNA-seq_and_meta_data/share_files/all_cells_count_matrix_filtered.tsv"
      )
    ),
    metadata = read_tsv(
      file.path(
        ag_data_dir,
        "Developmental_heart_filtered_scRNA-seq_and_meta_data/share_files/all_cells_meta_data_filtered.tsv"
      )
    )
  ),
  ST = list(
    count_matrix = read_tsv(
      file.path(
        ag_data_dir,
        "Developmental_heart_filtered_ST_matrix_and_meta_data/filtered_ST_matrix_and_meta_data/filtered_matrix.tsv"
      )
    ),
    metadata = read_tsv(
      file.path(
        ag_data_dir,
        "Developmental_heart_filtered_ST_matrix_and_meta_data/filtered_ST_matrix_and_meta_data/meta_data.tsv"
      )
    )
  )
)






co_meta_df <- readRDS(
  file.path(
    "data/scrna_heart_development",
    "descartes_Cao_ODay/df_cell.RDS"))

co_gene_df <- readRDS(
  file.path(
    "data/scrna_heart_development",
    "descartes_Cao_ODay/df_gene.RDS"))

co_count_matrix <- readRDS(
  file.path(
    "data/scrna_heart_development",
    "descartes_Cao_ODay/Heart_gene_count.RDS"))





stopifnot(identical(0L, sum(is.na(co_meta_df$sample))))
stopifnot(identical(
  co_meta_df$sample, unique(co_meta_df$sample)))


co_heart_meta_tbl <- as_tibble(co_meta_df) %>%
  filter(Organ == "Heart", !is.na(Main_cluster_name))

stopifnot(identical(0L, sum(is.na(co_heart_meta_tbl$Organ_cell_lineage))))

stopifnot(identical(0L, sum(is.na(co_heart_meta_tbl$sample))))
stopifnot(identical(
  co_heart_meta_tbl$sample, unique(co_heart_meta_tbl$sample)))

co_heart_meta_df <- column_to_rownames(co_heart_meta_tbl, "sample")

stopifnot(identical(0L, sum(is.na(rownames(co_heart_meta_df)))))
stopifnot(identical(
  rownames(co_heart_meta_df), unique(rownames(co_heart_meta_df))))
stopifnot(identical(
  co_heart_meta_tbl$sample, rownames(co_heart_meta_df)))


stopifnot(identical(0L, sum(is.na(co_gene_df))))
stopifnot(identical(
  rownames(co_gene_df), unique(rownames(co_gene_df))))

pp_co_gene_df <- co_gene_df
pp_co_gene_df$gene_id <- as.character(pp_co_gene_df$gene_id)
pp_co_gene_df$gene_type <- as.character(pp_co_gene_df$gene_type)
pp_co_gene_df$gene_short_name <- as.character(
  pp_co_gene_df$gene_short_name)

stopifnot(identical(0L, sum(is.na(pp_co_gene_df))))

stopifnot(identical(
  row.names(co_count_matrix),
  row.names(pp_co_gene_df)))
stopifnot(identical(
  pp_co_gene_df$gene_id,
  row.names(pp_co_gene_df)))



pp_co_count_matrix <- co_count_matrix[, row.names(co_heart_meta_df)]

stopifnot(identical(
  row.names(pp_co_count_matrix), row.names(pp_co_gene_df)))

stopifnot(identical(
  colnames(pp_co_count_matrix), row.names(co_heart_meta_df)))

pp_co_count_matrix <- pp_co_count_matrix[
  rowSums(pp_co_count_matrix) > 0, ]

stopifnot(identical(
  row.names(pp_co_count_matrix),
  unique(row.names(pp_co_count_matrix))))

f_pp_co_gene_df <- pp_co_gene_df[row.names(pp_co_count_matrix), ]
stopifnot(identical(
  row.names(f_pp_co_gene_df), row.names(pp_co_count_matrix)))
stopifnot(identical(
  row.names(f_pp_co_gene_df), f_pp_co_gene_df$gene_id))

row.names(pp_co_count_matrix) <- make.unique(
  f_pp_co_gene_df$gene_short_name)



stopifnot(identical(colnames(pp_co_count_matrix), row.names(co_heart_meta_df)))
co_seu <- CreateSeuratObject(
  pp_co_count_matrix,
  meta.data = co_heart_meta_df)

# The preprocessing procedure of the co dataset is the same as
# the one in
# <https://github.com/satijalab/azimuth-references/blob/b8b07dcdfcc09816a85aad07362e7bad4de03976/human_fetus/scripts/preprocess_BBI_reference.R>

colSums(is.na(co_seu@meta.data))

stopifnot(identical(
  as.character(co_seu@meta.data$Organ_cell_lineage),
  paste(
    as.character(co_seu@meta.data$Organ),
    as.character(co_seu@meta.data$Main_cluster_name),
    sep = "-")))

# According to <https://github.com/satijalab/azimuth-references/blob/b8b07dcdfcc09816a85aad07362e7bad4de03976/human_fetus/scripts/make_BBI_reference.R#L18-L30>:
# - l1 is Organ_cell_lineage without Organ prefix
# - l2 is Organ_cell_lineage
# - organ is the organ prefix in Organ_cell_lineage
#
# Validated in azimuth_fetus_ref.
azimuth_fetus_ref_sd <- SeuratData::LoadData(
  "fetusref", type = "azimuth")

azimuth_fetus_ref <- azimuth_fetus_ref_sd$map
stopifnot(identical(
  as.character(azimuth_fetus_ref@meta.data$annotation.l2),
  paste(
    as.character(azimuth_fetus_ref@meta.data$organ),
    as.character(azimuth_fetus_ref@meta.data$annotation.l1),
    sep = "-")))

pp_seu_list <- list(
  kd = CreateSeuratObject(
    kd_rds@assays$RNA$counts,
    meta.data = kd_rds@meta.data),

  cz = CreateSeuratObject(
    column_to_rownames(cz_tbl, "Gene")),

  ag_scRNA = CreateSeuratObject(
    column_to_rownames(
      ag_list$scRNA$count_matrix, "...1")),

  co = co_seu
)


saveRDS(pp_seu_list, file.path(out_dir, "pp_seu_list.rds"))
