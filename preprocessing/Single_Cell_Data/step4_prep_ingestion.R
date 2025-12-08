library(Seurat)
library(tidyverse)



analysis_id <- paste0(
  "prepare_ingestion_v1"
)
out_dir <- paste0("results/", analysis_id, "_results")
out_int_dir <- paste0("results/", analysis_id, "_intermediate_files")

if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

if (!dir.exists(out_int_dir)) {
  dir.create(out_int_dir)
}


co_seu <- readRDS(
  file.path(
    "results/prepare_azimuth_reference_v1_intermediate_files",
    "intermediate_bbi.rds"))

az_list <- readRDS(
  file.path(
    "results/azimuth_v2_intermediate_files",
    "az_list.rds"))



#- co
get_co_seu <- function() {
  rc_co_seu <- co_seu
  stopifnot(identical(
    0L,
    sum(is.na(
      rc_co_seu@meta.data[, c("Fetus_id",
                              "Development_day",
                              "Main_cluster_name",
                              "Organ_cell_lineage")]))
  ))
  rc_co_seu$Donor_id <- as.character(rc_co_seu$Fetus_id)
  rc_co_seu$Sample_id <- as.character(rc_co_seu$Fetus_id)
  rc_co_seu$Protocol <- "sci-RNA-seq3"
  rc_co_seu$Protocol_id <- "sci-RNA-seq3"
  rc_co_seu$Development_week <- paste0(
    round(rc_co_seu$Development_day / 7, 3),
    " Post-Conception Week")
  rc_co_seu$Development_week_short <- paste0(
    round(rc_co_seu$Development_day / 7, 3),
    "-PCW")
  rc_co_seu$Development_stage <- "fetal"
  rc_co_seu$Development_id <- paste(
    rc_co_seu$Development_stage,
    rc_co_seu$Development_week_short,
    sep = "-")

  rc_co_seu$Anatomical_annotation <- ""
  rc_co_seu$Azimuth_human_fetal_heart_dev_annotation.l1 <- as.character(
    rc_co_seu$Main_cluster_name)
  rc_co_seu$Azimuth_human_fetal_heart_dev_annotation.l2 <- as.character(
    rc_co_seu$Organ_cell_lineage)

  DefaultAssay(rc_co_seu) <- "RNA"
  rc_co_seu@reductions$ref.umap <- rc_co_seu@reductions$umap
  rc_co_seu@reductions$integrated_dr <- rc_co_seu@reductions$pca

  stopifnot(identical(
    rc_co_seu@reductions$ref.umap, rc_co_seu@reductions$umap))
  stopifnot(identical(
    rc_co_seu@reductions$integrated_dr, rc_co_seu@reductions$pca))

  return(rc_co_seu)
}



metadata_cols <- c(
  "Donor_id", "Sample_id",
  "Protocol", "Protocol_id",
  "Development_week", "Development_week_short",
  "Development_stage", "Development_id",
  "Anatomical_annotation", "Organ",
  "Azimuth_human_fetal_heart_dev_annotation.l1",
  "Azimuth_human_fetal_heart_dev_annotation.l2")


#- kd
get_kd_seu <- function() {
  rc_kd_seu <- az_list$kd
  stopifnot(identical(
    0L,
    sum(is.na(
      rc_kd_seu@meta.data[, c("sample", "donor", "region",
                              "stage",
                              "predicted.annotation.l1",
                              "predicted.annotation.l2")]))
  ))

  rc_kd_seu$Donor_id <- as.character(rc_kd_seu$donor)
  rc_kd_seu$Sample_id <- as.character(rc_kd_seu$sample)
  rc_kd_seu$Protocol <- "10x Chromium Single Cell 3' v3"
  rc_kd_seu$Protocol_id <- "10x-3'-v3"


  rc_kd_pcw <- round(rc_kd_seu$age.daysPC / 7, 3)
  stopifnot(identical(
    is.na(rc_kd_pcw),
    is.na(rc_kd_seu$age.daysPC)))


  rc_kd_seu$Development_week <- unname(if_else(
    is.na(rc_kd_pcw),
    true = "",
    false = paste0(rc_kd_pcw, " Post-Conception Week")))

  stopifnot(identical(
    rc_kd_seu$Development_week == "",
    is.na(rc_kd_pcw)))


  rc_kd_seu$Development_week_short <- unname(if_else(
    is.na(rc_kd_pcw),
    true = "",
    false = paste0(rc_kd_pcw, "-PCW")))

  stopifnot(identical(
    rc_kd_seu$Development_week_short == "",
    rc_kd_seu$Development_week == ""
  ))


  rc_kd_seu$Development_stage <- unname(if_else(
    rc_kd_seu$stage == "foetal",
    true = "fetal",
    false = "adult"))


  dev_tbl <- as_tibble(
    rc_kd_seu@meta.data, rownames = "cell_barcode"
  ) %>%
    select(
      cell_barcode, Development_week_short,
      Development_stage, age_group) %>%
    replace_na(list(age_group = ""))

  stopifnot(identical(0L, sum(is.na(dev_tbl))))
  stopifnot(identical(
    dev_tbl$Development_week_short == "",
    dev_tbl$Development_stage == "adult"))
  stopifnot(identical(
    !(dev_tbl$age_group == ""),
    dev_tbl$Development_stage == "adult"))

  dev_tbl <- dev_tbl %>%
    mutate(
      Development_id = if_else(
        dev_tbl$Development_stage == "adult",
        true = paste(
          dev_tbl$Development_stage, dev_tbl$age_group,
          "YO", sep = "-"),
        false = paste(
          dev_tbl$Development_stage,
          dev_tbl$Development_week_short,
          sep = "-")))

  stopifnot(identical(0L, sum(is.na(dev_tbl))))
  stopifnot(identical(
    pull(dev_tbl, cell_barcode),
    row.names(rc_kd_seu@meta.data)))

  rc_kd_seu$Development_id <- pull(
    dev_tbl, Development_id, name = cell_barcode)

  rc_kd_seu$Anatomical_annotation <- rc_kd_seu$region
  rc_kd_seu$Organ <- "Heart"
  rc_kd_seu$Azimuth_human_fetal_heart_dev_annotation.l1 <- as.character(
    rc_kd_seu$predicted.annotation.l1)
  rc_kd_seu$Azimuth_human_fetal_heart_dev_annotation.l2 <- as.character(
    rc_kd_seu$predicted.annotation.l2)

  return(rc_kd_seu)
}



#- cz
get_cz_seu <- function() {
  cz_barcode_tbl <- read_tsv(
    file.path(
      "data/scrna_heart_development/GSE106118_Cui_Zheng",
      "GSE106118_barcode_information.txt"
    )
  ) %>%
    mutate(
      cell_pfx = str_match(
        cell, "^([^.]+)\\.[0-9]+$")[, 2],
      file_pfx = str_match(
        prcessed_file, "^([^.]+)[_.]TPM\\.txt$")[, 2]) %>%
    select(!c(barcode))

  stopifnot(identical(0L, sum(is.na(cz_barcode_tbl))))


  cz_cell_file_pfx_tbl <- cz_barcode_tbl %>%
    count(file_pfx, cell_pfx) %>%
    arrange(file_pfx)





  rc_cz_seu <- az_list$cz
  stopifnot(identical(
    0L,
    sum(is.na(
      rc_cz_seu@meta.data[, c("predicted.annotation.l1",
                              "predicted.annotation.l2")]))
  ))
  stopifnot(identical(0L, sum(is.na(row.names(rc_cz_seu@meta.data)))))

  seu_cells_no_info <- row.names(rc_cz_seu@meta.data)[
    !(row.names(rc_cz_seu@meta.data) %in% cz_barcode_tbl$cell)]

  seu_cells_no_info_prev <- paste0(
    substr(seu_cells_no_info, 1, 10),
    as.character(as.numeric(substr(seu_cells_no_info, 11, 12)) - 1))

  seu_cells_no_info_prev_tbl <- cz_barcode_tbl %>%
    filter(cell %in% seu_cells_no_info_prev)

  stopifnot(identical(
    pull(seu_cells_no_info_prev_tbl, cell_pfx),
    pull(seu_cells_no_info_prev_tbl, file_pfx)))

  c_cz_barcode_tbl <- bind_rows(
    select(cz_barcode_tbl, prcessed_file, cell),
    tibble(
      prcessed_file = paste0(
        str_match(
          seu_cells_no_info, "^([^.]+)\\.[0-9]+$")[, 2],
        ".TPM.txt"),
      cell = seu_cells_no_info
    )
  ) %>%
    mutate(Sample_id = str_match(
      prcessed_file, "^([^.]+)[_.]TPM\\.txt$")[, 2])

  stopifnot(identical(0L, sum(is.na(c_cz_barcode_tbl))))

  stopifnot(identical(
    c_cz_barcode_tbl %>%
      distinct(prcessed_file, Sample_id) %>%
      pull(Sample_id),
    c_cz_barcode_tbl %>%
      distinct(Sample_id) %>%
      pull(Sample_id)
  ))

  stopifnot(identical(
    pull(c_cz_barcode_tbl, cell),
    unique(pull(c_cz_barcode_tbl, cell))))

  cz_cell_sample_lookup <- pull(
    c_cz_barcode_tbl, Sample_id, name = cell)

  stopifnot(identical(0L, sum(is.na(cz_cell_sample_lookup))))
  stopifnot(identical(
    pull(c_cz_barcode_tbl, cell),
    names(cz_cell_sample_lookup)))



  rc_cz_seu$cell_pfx <- str_match(
    row.names(rc_cz_seu@meta.data),
    "^([^.]+)\\.[0-9]+$")[, 2]

  rc_cz_seu$Anatomical_annotation <- str_match(
    row.names(rc_cz_seu@meta.data),
    "^HE[0-9]+W_[0-9]+_([A-Z]+)\\.[0-9]+$")[, 2]


  stopifnot(all(
    row.names(rc_cz_seu@meta.data) %in%
      names(cz_cell_sample_lookup)
  ))

  stopifnot(identical(
    row.names(rc_cz_seu@meta.data),
    unique(row.names(rc_cz_seu@meta.data))))

  rc_cz_seu$Donor_id <- str_match(
    row.names(rc_cz_seu@meta.data),
    "^(HE[0-9]+W_[0-9]+)_")[, 2]

  rc_cz_seu$Sample_id <- cz_cell_sample_lookup[
    row.names(rc_cz_seu@meta.data)]

  rc_cz_seu$Protocol <- "STRT-seq"
  rc_cz_seu$Protocol_id <- "STRT-seq"

  dev_week_match <- str_match(
    row.names(rc_cz_seu@meta.data),
    "^HE([0-9]+)W_[0-9]+_")[, 2]

  rc_cz_seu$Development_week <- paste0(
    dev_week_match, " Gestational Week")

  rc_cz_seu$Development_week_short <- paste0(
    dev_week_match, "-GW")

  rc_cz_seu$Development_stage <- "fetal"

  rc_cz_seu$Development_id <- paste(
    rc_cz_seu$Development_stage,
    rc_cz_seu$Development_week_short,
    sep = "-")

  rc_cz_seu$Anatomical_annotation <- str_match(
    row.names(rc_cz_seu@meta.data),
    "^HE[0-9]+W_[0-9]+_([^.]+)\\.[0-9]+$")[, 2]


  rc_cz_seu$Organ <- "Heart"
  rc_cz_seu$Azimuth_human_fetal_heart_dev_annotation.l1 <- as.character(
    rc_cz_seu$predicted.annotation.l1)
  rc_cz_seu$Azimuth_human_fetal_heart_dev_annotation.l2 <- as.character(
    rc_cz_seu$predicted.annotation.l2)

  stopifnot(identical(0L, sum(is.na(rc_cz_seu@meta.data))))

  return(rc_cz_seu)
}





#- ag
get_ag_seu <- function() {
  ag_meta_lines <- read_lines(
    file.path(
      "data/scrna_heart_development",
      "EGAS00001003996_asp_giacomello",
      "dev_heart_count_matrices_and_meta/Filtered",
      "Developmental_heart_filtered_scRNA-seq_and_meta_data",
      "share_files/all_cells_meta_data_filtered.tsv"))

  ag_meta_match <- str_match(
    ag_meta_lines,
    "^([ACGT]+\\.[0-9])\t[0-9]+\t[0-9]+\t(Exp_[0-9])\t[^\\t]+\t([0-9]+)\t(.+)$")

  ag_cell_exp_matrix <- ag_meta_match[!is.na(ag_meta_match[, 1]), 2:5]
  colnames(ag_cell_exp_matrix) <- c("cell", "Sample_id", "cluster", "end")

  ag_cell_exp_tbl <- as_tibble(ag_cell_exp_matrix)
  stopifnot(identical(0L, sum(is.na(ag_cell_exp_tbl))))
  stopifnot(identical(
    pull(ag_cell_exp_tbl, cell),
    unique(pull(ag_cell_exp_tbl, cell))))

  ag_cell_exp_lookup <- pull(ag_cell_exp_tbl, Sample_id, name = cell)
  stopifnot(identical(
    names(ag_cell_exp_lookup), pull(ag_cell_exp_tbl, cell)))


  rc_ag_seu <- az_list$ag_scRNA

  stopifnot(identical(
    0L,
    sum(is.na(
      rc_ag_seu@meta.data[, c("predicted.annotation.l1",
                              "predicted.annotation.l2")]))
  ))
  stopifnot(identical(0L, sum(is.na(row.names(rc_ag_seu@meta.data)))))

  stopifnot(all(
    row.names(rc_ag_seu@meta.data) %in%
      names(ag_cell_exp_lookup)))

  rc_ag_seu$Donor_id <- "Fetus1"
  rc_ag_seu$Sample_id <- ag_cell_exp_lookup[
    row.names(rc_ag_seu@meta.data)]
  rc_ag_seu$Protocol <- "10x Chromium Single Cell 3' v2"
  rc_ag_seu$Protocol_id <- "10x-3'-v2"

  rc_ag_seu$Development_week <- "6.5-7 Post-Conception Week"
  rc_ag_seu$Development_week_short <- "6.5-7-PCW"
  rc_ag_seu$Development_stage <- "fetal"
  rc_ag_seu$Development_id <- paste(
    rc_ag_seu$Development_stage,
    rc_ag_seu$Development_week_short,
    sep = "-")

  rc_ag_seu$Anatomical_annotation <- unname(if_else(
    rc_ag_seu$Sample_id == "Exp_1",
    true = "fraction (i)", false = "fraction (ii)"))
  rc_ag_seu$Organ <- "Heart"
  rc_ag_seu$Azimuth_human_fetal_heart_dev_annotation.l1 <- as.character(
    rc_ag_seu$predicted.annotation.l1)
  rc_ag_seu$Azimuth_human_fetal_heart_dev_annotation.l2 <- as.character(
    rc_ag_seu$predicted.annotation.l2)


  stopifnot(identical(0L, sum(is.na(rc_ag_seu@meta.data))))
  return(rc_ag_seu)
}


seu_list <- list(
  kd = get_kd_seu(),
  cz = get_cz_seu(),
  ag = get_ag_seu(),
  co = get_co_seu()
)


fmt_seu_list <- imap(
  list(
    co = "Cao_2020",
    kd = "Knight-Schrijver_2022",
    cz = "Cui_2019",
    ag = "Asp_2019"
  ),
  function(ds_id, xname) {
    print(xname)

    xseu <- seu_list[[xname]]
    stopifnot(!is.null(xseu))

    stopifnot(identical(0L, sum(is.na(xseu@meta.data[, metadata_cols]))))

    xseu$Dataset_id <- ds_id

    stopifnot(identical(0L, sum(is.na(xseu@meta.data[, metadata_cols]))))

    obj <- RenameCells(
      xseu, add.cell.id = paste0(ds_id, "__"))

    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj)
    obj <- ScaleData(obj)

    obj[["percent.mt"]] <- PercentageFeatureSet(
      obj, pattern = "^MT-")

    # Ref:
    # https://github.com/samuel-marsh/scCustomize/
    #   blob/e7efdaf18f2089b855569e8ff9796d299aa5c45a/R/Object_Utilities.R#L47-L190
    obj[["percent.ribo"]] <- PercentageFeatureSet(
      obj, pattern = "^RP[SL]")

    obj[["percent.hb"]] <- PercentageFeatureSet(
      obj, pattern = "^HB[^(P)]")

    obj <- subset(obj, subset = nCount_RNA > 250)

    cat("--------------------------\n\n\n")
    return(obj)
  }
)

print("nCount_RNA quantiles")
print(map(seu_list, function(x) quantile(x$nCount_RNA, seq(0, 0.25, 0.05))))
print("Cell counts")
print(map(seu_list, function(x) table(x$Development_stage)))
print("Number of genes with > 0 counts in > 0 cells")
print(map(seu_list, function(x) sum(rowSums(x@assays$RNA$counts > 0) > 0)))


print("After filtering")
print("nCount_RNA quantiles")
print(map(fmt_seu_list, function(x) quantile(x$nCount_RNA, seq(0, 0.25, 0.05))))
print("Cell counts")
print(map(fmt_seu_list, function(x) table(x$Development_stage)))
print("Number of genes with > 0 counts in > 0 cells")
print(map(fmt_seu_list, function(x) sum(rowSums(x@assays$RNA$counts > 0) > 0)))

save_file <- file.path(out_dir, "fmt_seu_list.rds")
saveRDS(fmt_seu_list, save_file)

save_file_sha256 <- tools::sha256sum(save_file)
save_file_sha256_line <- paste(
  save_file_sha256, names(save_file_sha256),
  sep = "  ")

write_lines(
  save_file_sha256_line,
  file.path(out_dir, "sha256sum.txt"))
