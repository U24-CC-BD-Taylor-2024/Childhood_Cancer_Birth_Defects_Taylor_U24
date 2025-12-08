library(Seurat)
library(tidyverse)

RESERVED_CHAR_REPLACEMENT <- "."
SC_RNA_CLUST_SAB <- "HDSCCEXP"

stopifnot(identical(
  length(RESERVED_CHAR_REPLACEMENT), 1L))
stopifnot(!is.na(RESERVED_CHAR_REPLACEMENT))
stopifnot(identical(
  length(SC_RNA_CLUST_SAB), 1L))
stopifnot(!is.na(SC_RNA_CLUST_SAB))


analysis_id <- paste0(
  "seurat_to_owlnets_v1"
)
out_dir <- paste0("results/", analysis_id, "_results")
out_int_dir <- paste0(
  "results/", analysis_id, "_intermediate_files")

if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

if (!dir.exists(out_int_dir)) {
  dir.create(out_int_dir)
}






fmt_seu_list <- readRDS(
  file.path(
    "results/prepare_ingestion_v1_results",
    "fmt_seu_list.rds"))

gsb_hgnc_tbl <- read_tsv(
  "data/hgnc_master.txt",
  guess_max = 1e6)

mapping_gsb_hgnc_tbl <- gsb_hgnc_tbl %>%
  select(symbol, hgnc_id) %>%
  rename(gene_symbol = symbol)


assert_one_to_one_mapping <- function(x, y) {
  stopifnot(is.vector(x))
  stopifnot(is.vector(y))
  stopifnot(identical(0L, sum(is.na(x))))
  stopifnot(identical(0L, sum(is.na(y))))
  stopifnot(identical(length(x), length(y)))

  tbl <- distinct(tibble(a = x, b = y))

  stopifnot(identical(pull(tbl, a), unique(pull(tbl, a))))
  stopifnot(identical(pull(tbl, b), unique(pull(tbl, b))))
}

assert_many_to_one_mapping <- function(x, y) {
  stopifnot(is.vector(x))
  stopifnot(is.vector(y))
  stopifnot(identical(0L, sum(is.na(x))))
  stopifnot(identical(0L, sum(is.na(y))))
  stopifnot(identical(length(x), length(y)))

  tbl <- distinct(tibble(a = x, b = y))

  stopifnot(identical(pull(tbl, a), unique(pull(tbl, a))))
}

stopifnot(identical(0L, sum(is.na(mapping_gsb_hgnc_tbl))))
stopifnot(identical(
  mapping_gsb_hgnc_tbl,
  distinct(mapping_gsb_hgnc_tbl)
))
assert_one_to_one_mapping(
  pull(mapping_gsb_hgnc_tbl, gene_symbol),
  pull(mapping_gsb_hgnc_tbl, hgnc_id)
)




RUN_DATASETS <- c("co", "kd", "cz", "ag")

run_seu_list <- map(fmt_seu_list[RUN_DATASETS], function(xseu) {
  count_matrix <- GetAssayData(
    GetAssay(xseu, assay = "RNA"),
    layer = "counts")

  stopifnot(identical(row.names(count_matrix), rownames(xseu)))
  stopifnot(identical(colnames(count_matrix), colnames(xseu)))

  filter_genes <- row.names(count_matrix)[
    rowSums(count_matrix > 0) > 0]

  res <- subset(xseu, features = filter_genes)

  return(res)
})



annot_stats_tbl <- bind_rows(
  map(run_seu_list, function(xseu) {
    meta_tbl <- xseu@meta.data %>%
      as_tibble() %>%
      count(
        Dataset_id,
        Azimuth_human_fetal_heart_dev_annotation.l2,
        Development_id,
        name = "n_cells") %>%
      mutate(
        percent_cells_in_dataset = (
          n_cells / sum(n_cells) * 100
        )
      )%>%
      mutate(
        dev_num = as.numeric(str_match(
          Development_id,
          "^[^-]+-([^-]+)-")[, 2]),
        dev_stage = str_match(
          Development_id,
          "^([^-]+)-[^-]+-")[, 2]
      )

    stopifnot(identical(0L, sum(is.na(meta_tbl))))

    res <- meta_tbl

    return(res)
  })
) %>%
  arrange(
    Azimuth_human_fetal_heart_dev_annotation.l2,
    Dataset_id,
    desc(dev_stage), dev_num
  ) %>%
  select(
    Azimuth_human_fetal_heart_dev_annotation.l2,
    Dataset_id, Development_id, n_cells,
    percent_cells_in_dataset
  )


write_tsv(
  annot_stats_tbl,
  file.path(
    out_dir,
    paste0(
      "scRNA_seq",
      "_manual_annotation_related_",
      "cellType_DevelopmentStage_stats.tsv")))




type_wide_list <- list(
  normalized_umi_count = map(run_seu_list, function(xseu) {
    wide_sparse <- GetAssayData(
      GetAssay(xseu, assay = "RNA"),
      layer = "data")

    # (cells, genes)
    wide_matrix <- t(expm1(as.matrix(wide_sparse)))
    stopifnot(identical(0L, sum(is.na(row.names(wide_matrix)))))
    stopifnot(identical(0L, sum(is.na(colnames(wide_matrix)))))
    stopifnot(identical(colnames(wide_matrix), row.names(xseu)))
    stopifnot(identical(row.names(wide_matrix), colnames(xseu)))
    stopifnot(identical(unique(colnames(wide_matrix)), row.names(xseu)))
    stopifnot(identical(unique(row.names(wide_matrix)), colnames(xseu)))
    stopifnot(all.equal(
      unname(rowSums(wide_matrix)),
      rep(10000, nrow(wide_matrix))))

    value_range <- c(min(wide_matrix), max(wide_matrix))

    wide_tbl <- as_tibble(wide_matrix, rownames = "cell_barcode")
    stopifnot(identical(pull(wide_tbl, cell_barcode), colnames(xseu)))
    stopifnot(identical(
      colnames(wide_tbl), c("cell_barcode", row.names(xseu))))
    stopifnot(identical(0L, sum(is.na(wide_tbl))))


    meta_tbl <- as_tibble(xseu@meta.data, rownames = "cell_barcode")
    stopifnot(identical(
      pull(wide_tbl, cell_barcode),
      pull(meta_tbl, cell_barcode)))

    res <- list(
      wide_tbl = wide_tbl,
      meta_tbl = meta_tbl,
      value_range = value_range
    )

    return(res)
  }),

  umap_embedding = map(run_seu_list, function(xseu) {
    wide_matrix0 <- Embeddings(xseu, "ref.umap")
    wide_matrix <- wide_matrix0
    row.names(wide_matrix) <- unname(row.names(wide_matrix))
    stopifnot(all(row.names(wide_matrix) == row.names(wide_matrix0)))
    stopifnot(identical(
      str_to_lower(colnames(wide_matrix)),
      c("umap_1", "umap_2")
    ))
    colnames(wide_matrix) <- c("Dim1", "Dim2")
    stopifnot(identical(row.names(wide_matrix), colnames(xseu)))
    stopifnot(identical(0L, sum(is.na(row.names(wide_matrix)))))
    stopifnot(identical(unique(row.names(wide_matrix)), colnames(xseu)))

    value_range <- c(min(wide_matrix), max(wide_matrix))

    wide_tbl <- as_tibble(wide_matrix, rownames = "cell_barcode")
    stopifnot(identical(
      colnames(wide_tbl),
      c("cell_barcode", "Dim1", "Dim2")
    ))
    stopifnot(identical(
      pull(wide_tbl, cell_barcode), colnames(xseu)))
    stopifnot(identical(0L, sum(is.na(wide_tbl))))

    meta_tbl <- as_tibble(xseu@meta.data, rownames = "cell_barcode")
    stopifnot(identical(
      pull(wide_tbl, cell_barcode),
      pull(meta_tbl, cell_barcode)))

    res <- list(
      wide_tbl = wide_tbl,
      meta_tbl = meta_tbl,
      value_range = value_range
    )

    return(res)
  })
)


type_range_list <- map(type_wide_list, function(one_type) {
  one_type_range_list <- map(one_type, function(wlist) {
    res <- wlist$value_range

    return(res)
  })

  reduce(one_type_range_list, function(x, y) {
    stopifnot(is.numeric(x))
    stopifnot(is.numeric(y))
    stopifnot(identical(length(x), 2L))
    stopifnot(identical(length(y), 2L))
    stopifnot(identical(0L, sum(is.na(x))))
    stopifnot(identical(0L, sum(is.na(y))))


    res <- c(min(c(x[1], y[1])), max(c(x[2], y[2])))

    stopifnot(is.numeric(res))
    stopifnot(identical(length(res), 2L))
    return(res)
  })
})


type_sab_lookup <- c(
  normalized_umi_count = "SCNORMBINS",
  umap_embedding = "EMBREGIONS"
)

type_bin_config_list <- list(
  normalized_umi_count = list(
    lb = 0, ub = 10000, bw = 1, ndims = 1L
  ),
  umap_embedding = list(
    lb = floor(type_range_list$umap_embedding[1]),
    ub = ceiling(type_range_list$umap_embedding[2]),
    bw = 0.1, ndims = 2L
  )
)


bin_config_to_breaks <- function(xconfig) {
    stopifnot(is.numeric(xconfig$lb))
    stopifnot(is.numeric(xconfig$ub))
    stopifnot(is.numeric(xconfig$bw))
    stopifnot(is.integer(xconfig$ndims))

    stopifnot(!is.na(xconfig$lb))
    stopifnot(!is.na(xconfig$ub))
    stopifnot(!is.na(xconfig$bw))
    stopifnot(!is.na(xconfig$ndims))

    stopifnot(xconfig$lb < xconfig$ub)
    stopifnot(xconfig$bw > 0)
    stopifnot(xconfig$lb + xconfig$bw < xconfig$ub)
    stopifnot(xconfig$ndims > 0)
    stopifnot(identical(xconfig$bw, round(xconfig$bw, 1)))

    break_seq <- seq(from = xconfig$lb, to = xconfig$ub, by = xconfig$bw)

    stopifnot(identical(break_seq[length(break_seq)], xconfig$ub))
    stopifnot(identical(break_seq[1], xconfig$lb))

    # bins are [lb, ub), so add an additional bin to the right
    bin_breaks <- c(break_seq, xconfig$ub + xconfig$bw)

    round_bin_breaks <- round(bin_breaks, 3)
    stopifnot(all.equal(bin_breaks, round_bin_breaks))
    stopifnot(is.numeric(round_bin_breaks))

    stopifnot(identical(0L, sum(is.na(round_bin_breaks))))

    return(round_bin_breaks)
}

type_bin_breaks_list <- map(
  type_bin_config_list,
  bin_config_to_breaks
)


# conventional cut
cvn_cut <- function(x_vals, x_breaks) {
  stopifnot(is.numeric(x_vals))
  stopifnot(is.numeric(x_breaks))
  stopifnot(identical(0L, sum(is.na(x_vals))))
  stopifnot(identical(0L, sum(is.na(x_breaks))))
  stopifnot(identical(x_breaks, unique(x_breaks)))

  res <- cut(x_vals, breaks = x_breaks, right = FALSE)

  stopifnot(identical(0L, sum(is.na(res))))

  return(res)
}


type_odb_list <- imap(
  type_bin_config_list,
  function(xconfig, xtype) {
    x_breaks <- type_bin_breaks_list[[xtype]]
    stopifnot(is.numeric(x_breaks))

    x_one_dim_bins <- levels(cvn_cut(0, x_breaks))

    stopifnot(is.character(x_one_dim_bins))
    stopifnot(identical(0L, sum(is.na(x_one_dim_bins))))
    stopifnot(identical(
      x_one_dim_bins, unique(x_one_dim_bins)))

    stopifnot(identical(
      length(x_one_dim_bins),
      length(x_breaks) - 1L))


    one_dim_tbl <- tibble(code = x_one_dim_bins)

    stopifnot(all(!str_detect(pull(one_dim_tbl, code), "e")))
    stopifnot(identical(
      pull(one_dim_tbl, code),
      unique(pull(one_dim_tbl, code))
    ))

    one_dim_code_match <- str_match(
      pull(one_dim_tbl, code),
      "^\\[(-?[0-9]+\\.?[0-9]*),(-?[0-9]+\\.?[0-9]*)\\)$")

    stopifnot(identical(0L, sum(is.na(one_dim_code_match))))
    colnames(one_dim_code_match) <- c("code", "lb", "ub")

    check_tbl <- as_tibble(one_dim_code_match) %>%
      mutate(lb = as.numeric(lb), ub = as.numeric(ub)) %>%
      select(lb, ub, code)

    stopifnot(identical(0L, sum(is.na(check_tbl))))
    stopifnot(identical(
      pull(check_tbl, code),
      pull(one_dim_tbl, code)))

    write_tsv(
      check_tbl,
      file.path(
        out_int_dir,
        paste0(type_sab_lookup[xtype], "_one_dim_tbl.tsv")))

    return(check_tbl)
  }
)

type_link_out_node_list <- imap(
  type_bin_config_list,
  function(xconfig, xtype) {
    one_dim_tbl <- type_odb_list[[xtype]]
    stopifnot(is_tibble(one_dim_tbl))
    stopifnot(identical(0L, sum(is.na(one_dim_tbl))))

    if (xconfig$ndims > 1) {
      dim_inds <- seq(from = 1L, to = xconfig$ndims, by = 1L)
      stopifnot(identical(length(dim_inds), xconfig$ndims))
      dim_ps <- paste0("Dim", dim_inds)

      ndim_code_list <- map(set_names(dim_ps), function(x_dim) {
        res <- paste0(x_dim, pull(one_dim_tbl, code))
        stopifnot(identical(0L, sum(is.na(res))))
        stopifnot(is.character(res))
        return(res)
      })

      ndim_code_grid_tbl <- do.call(
        expand_grid, ndim_code_list)

      stopifnot(identical(0L, sum(is.na(ndim_code_grid_tbl))))

      ndim_code_tbl <- ndim_code_grid_tbl %>%
        unite("code", everything(), sep = ";")

    } else {
      stopifnot(identical(xconfig$ndims, 1L))

      ndim_code_tbl <- one_dim_tbl
    }

    ndim_code_tbl <- ndim_code_tbl %>%
      mutate(SAB = type_sab_lookup[xtype])

    stopifnot(identical(0L, sum(is.na(ndim_code_tbl))))
    stopifnot(identical(
      pull(ndim_code_tbl, code),
      unique(pull(ndim_code_tbl, code))
    ))

    res <- ndim_code_tbl %>%
      mutate(
        node_id = paste(SAB, code, sep = ":"),
        node_label = "",
        node_definition = "",
        node_synonyms = "",
        node_dbxrefs = "",
        node_namespace = ""
      ) %>%
      select(
        node_id, node_label, node_definition,
        node_synonyms, node_dbxrefs, node_namespace
      )

    stopifnot(identical(0L, sum(is.na(res))))
    stopifnot(identical(
      pull(res, node_id),
      unique(pull(res, node_id))
    ))

    return(res)
  }
)



iwalk(type_link_out_node_list, function(xtbl, xtype) {
  sab_out_dir <- file.path(out_dir, type_sab_lookup[xtype])

  if (!dir.exists(sab_out_dir)) {
    dir.create(sab_out_dir)
  }

  write_tsv(
    xtbl,
    file.path(sab_out_dir, "OWLNETS_node_metadata.txt"))

  empty_edge_tbl <- tibble(
    subject = character(0),
    predicate = character(0),
    object = character(0)
  )

  write_tsv(
    empty_edge_tbl,
    file.path(sab_out_dir, "OWLNETS_edgelist.txt"))
})






metadata_cols <- c(
  "Dataset_id", "Donor_id", "Sample_id",
  "Protocol", "Protocol_id",
  "Development_week", "Development_week_short",
  "Development_stage", "Development_id",
  "Anatomical_annotation", "Organ",
  "Azimuth_human_fetal_heart_dev_annotation.l1",
  "Azimuth_human_fetal_heart_dev_annotation.l2")

type_pb_wide_list <- imap(type_wide_list, function(xtype, xn) {
  imap(xtype, function(x_wide, wn) {
    print(xn)
    print(wn)
    stopifnot(identical(
      0L, sum(is.na(x_wide$meta_tbl[, c("cell_barcode", metadata_cols)]))
    ))
    stopifnot(identical(
      pull(x_wide$wide_tbl, cell_barcode),
      pull(x_wide$meta_tbl, cell_barcode)))



    cb_pb_tbl <- x_wide$meta_tbl %>%
      mutate(
        pb_group = paste(
          Dataset_id,
          Sample_id,
          Development_id,
          Protocol_id,
          Azimuth_human_fetal_heart_dev_annotation.l2,
          sep = "..."
        )
      ) %>%
      mutate(
        cb_group = str_replace_all(
          pb_group, "[ :_]", RESERVED_CHAR_REPLACEMENT
        )
      ) %>%
      select(cell_barcode, cb_group, pb_group)

    assert_one_to_one_mapping(
      pull(cb_pb_tbl, cb_group),
      pull(cb_pb_tbl, pb_group)
    )
    stopifnot(identical(0L, sum(is.na(pull(cb_pb_tbl, pb_group)))))
    stopifnot(is.character(pull(cb_pb_tbl, pb_group)))

    stopifnot(identical(0L, sum(is.na(pull(cb_pb_tbl, cb_group)))))
    stopifnot(is.character(pull(cb_pb_tbl, cb_group)))
    stopifnot(all(!str_detect(pull(cb_pb_tbl, cb_group), coll(" "))))
    stopifnot(all(!str_detect(pull(cb_pb_tbl, cb_group), coll(":"))))
    stopifnot(all(!str_detect(pull(cb_pb_tbl, cb_group), coll("_"))))

    cb_pb_tbl <- cb_pb_tbl %>%
      select(cell_barcode, cb_group)

    stopifnot(identical(0L, sum(is.na(cb_pb_tbl))))
    stopifnot(identical(
      cb_pb_tbl,
      distinct(cb_pb_tbl)))
    stopifnot(identical(
      pull(cb_pb_tbl, cell_barcode),
      unique(pull(cb_pb_tbl, cell_barcode))))



    stopifnot(identical(
      pull(x_wide$wide_tbl, cell_barcode),
      pull(cb_pb_tbl, cell_barcode)))

    cb_pb_wide_tbl <- x_wide$wide_tbl %>%
      mutate(
        cb_group = pull(cb_pb_tbl, cb_group)
      )

    stopifnot(identical(
      pull(cb_pb_wide_tbl, cell_barcode),
      pull(cb_pb_tbl, cell_barcode)))
    stopifnot(identical(
      pull(cb_pb_wide_tbl, cb_group),
      pull(cb_pb_tbl, cb_group)))




    pb_wide_tbl <- cb_pb_wide_tbl %>%
      select(!c(cell_barcode)) %>%
      group_by(cb_group) %>%
      group_modify(function(x_grp, x_key) {
        x_grp_matrix <- as.matrix(x_grp)
        stopifnot(identical(
          colnames(x_grp_matrix),
          colnames(x_grp)
        ))

        gs <- as_tibble_row(colMeans(x_grp_matrix)) %>%
          mutate(n_cells = nrow(x_grp_matrix)) %>%
          relocate(n_cells)

        stopifnot(identical(
          colnames(gs),
          c("n_cells", colnames(x_grp))
        ))

        stopifnot(identical(0L, sum(is.na(gs))))

        return(gs)
      }) %>%
      ungroup()

    print(quantile(pull(pb_wide_tbl, n_cells), seq(0, 1, 0.1)))

    return(pb_wide_tbl)
  })
})






ds_vn_id_lookup <- c(
  co = "Cao_2020",
  kd = "Knight-Schrijver_2022",
  cz = "Cui_2019",
  ag = "Asp_2019"
)



cb_group_to_tbl <- function(cb_group_vec) {
  id_separator <- "\\.\\.\\."

  cb_group_match <- str_match(
    cb_group_vec,
    paste0(
      "^(.+)",
      id_separator,
      "(.+)",
      id_separator,
      "(.+)",
      id_separator,
      "(.+)",
      id_separator,
      "(.+)$"
    )
  )
  stopifnot(identical(0L, sum(is.na(cb_group_match))))
  colnames(cb_group_match) <- c(
    "cb_group",
    "Dataset_id", "Sample_id", "Development_id",
    "Protocol_id", "Cluster_id")

  res <- as_tibble(cb_group_match) %>%
    mutate(
      Development_stage = str_match(
        Development_id, "^([^-]+)-")[, 2])

  stopifnot(identical(0L, sum(is.na(res))))
  return(res)
}


pb_stats_type <- "normalized_umi_count"

pb_stats_tbl <- bind_rows(
  imap(type_pb_wide_list[[pb_stats_type]], function(xtbl, xname) {
    stopifnot(identical(0L, sum(is.na(xtbl))))
    stopifnot(identical(
      pull(xtbl, cb_group),
      unique(pull(xtbl, cb_group))
    ))

    sd_tbl <- type_wide_list[[pb_stats_type]][[xname]]$meta_tbl %>%
      select(Sample_id, Donor_id) %>%
      distinct()

    assert_many_to_one_mapping(
      pull(sd_tbl, Sample_id),
      pull(sd_tbl, Donor_id)
    )

    sample_donor_tbl <- sd_tbl %>%
      mutate(Sample_id = str_replace_all(
        Sample_id, "[ :_]", RESERVED_CHAR_REPLACEMENT))

    stopifnot(identical(0L, sum(is.na(sample_donor_tbl))))
    stopifnot(identical(
      distinct(sample_donor_tbl, Sample_id),
      select(sample_donor_tbl, Sample_id)
    ))

    assert_one_to_one_mapping(
      pull(sd_tbl, Sample_id),
      pull(sample_donor_tbl, Sample_id)
    )


    cbg_tbl <- cb_group_to_tbl(pull(xtbl, cb_group))

    cbg_n_tbl <- cbg_tbl %>%
      left_join(
        select(xtbl, cb_group, n_cells),
        by = "cb_group",
        relationship = "one-to-one") %>%
      left_join(
        sample_donor_tbl,
        by = "Sample_id",
        relationship = "many-to-one")

    stopifnot(identical(cbg_n_tbl, distinct(cbg_n_tbl)))
    stopifnot(identical(0L, sum(is.na(cbg_n_tbl))))

    res <- cbg_n_tbl %>%
      group_by(Development_stage) %>%
      summarize(
        n_fetuses_or_donors = length(unique(Donor_id)),
        n_samples = length(unique(Sample_id)),
        n_clusters = length(unique(Cluster_id)),
        n_cells = sum(n_cells),
        n_pseudobulks = n()
      ) %>%
      ungroup() %>%
      mutate(
        n_genes = ncol(xtbl) - 2,
        dataset = ds_vn_id_lookup[xname]) %>%
      relocate(dataset) %>%
      arrange(dataset, Development_stage)

    stopifnot(nrow(res) %in% c(1L, 2L))
    stopifnot(identical(0L, sum(is.na(res))))

    return(res)
  })
)

write_tsv(
  pb_stats_tbl,
  file.path(
    out_dir,
    paste0(
      "scRNA_seq",
      "_data_source_stats.tsv")))

print(data.frame(pb_stats_tbl))
iwalk(type_wide_list[[pb_stats_type]], function(x, xname) {
  print(xname)
  print(length(unique(x$meta_tbl$Donor_id)))
  print(sort(unique(x$meta_tbl$Donor_id)))
})







type_pb_long_list <- map(type_pb_wide_list, function(xtype) {
  map(xtype, function(pb_wide_tbl) {
    x_long_tbl <- pb_wide_tbl %>%
      select(!c(n_cells)) %>%
      pivot_longer(
        !cb_group,
        names_to = "FEATURE_NAME",
        values_to = "FEATURE_VALUE")

    stopifnot(identical(0L, sum(is.na(x_long_tbl))))
    stopifnot(identical(
      colnames(x_long_tbl),
      c("cb_group", "FEATURE_NAME", "FEATURE_VALUE")))

    stopifnot(identical(
      sort(unique(pull(x_long_tbl, cb_group))),
      sort(pull(pb_wide_tbl, cb_group))
    ))

    return(x_long_tbl)
  })
})






validate_val_bin_mappings <- function(vals, bins, all_bins) {
  stopifnot(is.character(bins))
  stopifnot(is.character(all_bins))
  stopifnot(is.numeric(vals))
  stopifnot(identical(0L, sum(is.na(vals))))
  stopifnot(identical(0L, sum(is.na(bins))))
  stopifnot(identical(0L, sum(is.na(all_bins))))
  stopifnot(identical(all_bins, unique(all_bins)))
  stopifnot(identical(length(vals), length(bins)))

  bin_match <- str_match(
    bins,
    "^\\[(-?[0-9]+\\.?[0-9]*),(-?[0-9]+\\.?[0-9]*)\\)$")

  stopifnot(identical(0L, sum(is.na(bin_match))))
  lb <- as.numeric(bin_match[, 2])
  ub <- as.numeric(bin_match[, 3])

  stopifnot(identical(0L, sum(is.na(lb))))
  stopifnot(identical(0L, sum(is.na(ub))))

  stopifnot(all(vals >= lb))
  stopifnot(all(vals < ub))
  stopifnot(all(bins %in% all_bins))

  return(NULL)
}



primary_type <- "normalized_umi_count"

pt_long_list <- map(
  type_pb_long_list[[primary_type]],
  function(xtbl) {
    ptl_tbl <- xtbl %>%
      rename(
        gene_symbol = FEATURE_NAME,
        primary_val = FEATURE_VALUE
      ) %>%
      mutate(
        cb_group_code = paste0(
          SC_RNA_CLUST_SAB, ":", cb_group, "...",
          gene_symbol))


    stopifnot(identical(0L, sum(is.na(ptl_tbl))))
    stopifnot(identical(
      pull(ptl_tbl, cb_group_code),
      unique(pull(ptl_tbl, cb_group_code))))

    return(ptl_tbl)
  }
)


ds_node_tbl_list <- map(
  pt_long_list,
  function(xtbl) {
    node_tbl <- xtbl %>%
      select(cb_group_code) %>%
      rename(node_id = cb_group_code) %>%
      mutate(
        node_label = "",
        node_definition = "",
        node_synonyms = "",
        node_dbxrefs = "",
        node_namespace = ""
      ) %>%
      select(
        node_id, node_label, node_definition,
        node_synonyms, node_dbxrefs, node_namespace
      )

    stopifnot(identical(0L, sum(is.na(node_tbl))))
    stopifnot(identical(
      pull(node_tbl, node_id),
      unique(pull(node_tbl, node_id))))

    return(node_tbl)
  }
)



vals_to_bins <- function(val_vec, break_vec) {
  val_cut <- cvn_cut(val_vec, break_vec)

  val_bin_vec <- as.character(val_cut)
  stopifnot(identical(0L, sum(is.na(val_bin_vec))))
  stopifnot(all(val_bin_vec %in% levels(val_cut)))

  stopifnot(identical(
    length(val_bin_vec),
    length(val_vec)))
  stopifnot(identical(
    length(val_bin_vec),
    length(val_cut)))

  return(val_bin_vec)
}


ds_norm_uc_etbl_list <- imap(
  pt_long_list,
  function(x_long_tbl, xname) {
    edge_tbl <- x_long_tbl %>%
      rename(norm_uc = primary_val)

    stopifnot(identical(0L, sum(is.na(edge_tbl))))

    norm_uc_vec <- pull(edge_tbl, norm_uc)

    norm_uc_bin_vec <- vals_to_bins(
      norm_uc_vec,
      type_bin_breaks_list[["normalized_umi_count"]])


    edge_tbl$norm_uc_bin <- norm_uc_bin_vec
    stopifnot(identical(
      edge_tbl$norm_uc_bin, norm_uc_bin_vec))

    validate_val_bin_mappings(
      pull(edge_tbl, norm_uc),
      pull(edge_tbl, norm_uc_bin),
      pull(type_odb_list[["normalized_umi_count"]], code)
    )

    edge_tbl <- edge_tbl %>%
      mutate(
        norm_uc_code = paste0(
          type_sab_lookup["normalized_umi_count"],
          ":", norm_uc_bin))

    stopifnot(identical(0L, sum(is.na(edge_tbl))))

    edge_tbl <- edge_tbl %>%
      select(cb_group_code, norm_uc_code) %>%
      mutate(predicate = "has_expression") %>%
      rename(
        subject = cb_group_code,
        object = norm_uc_code) %>%
      select(subject, predicate, object)

    stopifnot(identical(0L, sum(is.na(edge_tbl))))
    stopifnot(all(
      pull(edge_tbl, object) %in%
        pull(
          type_link_out_node_list$normalized_umi_count,
          node_id)
    ))

    return(edge_tbl)
  }
)


ds_umap_etbl_list <- imap(
  type_pb_long_list[["umap_embedding"]],
  function(x_long_tbl, xname) {

    umap_tbl <- x_long_tbl %>%
      rename(
        dim_name = FEATURE_NAME,
        dim_val = FEATURE_VALUE) %>%
      mutate(
        dim_ind = as.numeric(
          str_match(dim_name, "^Dim([0-9]+)$")[, 2]))

    stopifnot(is_tibble(umap_tbl))
    stopifnot(identical(0L, sum(is.na(umap_tbl))))


    umap_vec <- pull(umap_tbl, dim_val)
    stopifnot(is.numeric(umap_vec))

    umap_bin_vec <- vals_to_bins(
      umap_vec,
      type_bin_breaks_list[["umap_embedding"]])

    umap_tbl$umap_bin <- umap_bin_vec

    validate_val_bin_mappings(
      pull(umap_tbl, dim_val),
      pull(umap_tbl, umap_bin),
      pull(type_odb_list[["umap_embedding"]], code)
    )
    stopifnot(identical(0L, sum(is.na(umap_tbl))))

    umap_cbg_annot_tbl <- umap_tbl %>%
      group_by(cb_group) %>%
      group_modify(function(x_grp, x_key) {
        gs <- x_grp %>%
          arrange(dim_ind) %>%
          mutate(dim_code = paste0(dim_name, umap_bin))

        gs_code <- paste(pull(gs, dim_code), collapse = ";")
        stopifnot(is.character(gs_code))
        stopifnot(identical(length(gs_code), 1L))

        res <- tibble(
          umap_code = paste0(
            type_sab_lookup["umap_embedding"],
            ":", gs_code))

        stopifnot(identical(nrow(res), 1L))

        return(res)
      }) %>%
      ungroup()

    stopifnot(identical(
      pull(umap_cbg_annot_tbl, cb_group),
      unique(pull(umap_cbg_annot_tbl, cb_group))))



    edge_tbl <- pt_long_list[[xname]] %>%
      left_join(
        umap_cbg_annot_tbl, by = "cb_group",
        relationship = "many-to-one")

    stopifnot(identical(
      pull(edge_tbl, cb_group),
      pull(pt_long_list[[xname]], cb_group)
    ))
    stopifnot(identical(
      pull(edge_tbl, cb_group_code),
      pull(pt_long_list[[xname]], cb_group_code)
    ))
    stopifnot(identical(0L, sum(is.na(edge_tbl))))

    edge_tbl <- edge_tbl %>%
      select(cb_group_code, umap_code) %>%
      mutate(predicate = "has_embedding") %>%
      rename(
        subject = cb_group_code,
        object = umap_code) %>%
      select(subject, predicate, object)

    stopifnot(all(
      pull(edge_tbl, object) %in%
        pull(
          type_link_out_node_list$umap_embedding,
          node_id)
    ))

    return(edge_tbl)
  }
)




ds_hgnc_annot_list <- imap(
  pt_long_list,
  function(x_long_tbl, xname) {
    hgnc_tbl <- x_long_tbl %>%
      left_join(
        mapping_gsb_hgnc_tbl %>%
          rename(hgnc_code = hgnc_id),
        by = "gene_symbol",
        relationship = "many-to-one")

    stopifnot(identical(
      pull(x_long_tbl, gene_symbol),
      pull(hgnc_tbl, gene_symbol)
    ))
    stopifnot(identical(
      pull(x_long_tbl, cb_group_code),
      pull(hgnc_tbl, cb_group_code)
    ))
    stopifnot(identical(
      0L,
      sum(is.na(select(hgnc_tbl, !c(hgnc_code))))
    ))


    return(hgnc_tbl)
  }
)

ds_gsb_no_hgnc_tbl_list <- imap(
  ds_hgnc_annot_list,
  function(hgnc_tbl, xname) {
    res <- hgnc_tbl %>%
      filter(is.na(hgnc_code)) %>%
      select(gene_symbol) %>%
      distinct() %>%
      arrange(gene_symbol)

    stopifnot(identical(0L, sum(is.na(res))))

    write_tsv(
      res,
      file.path(
        out_int_dir,
        paste0(
          xname,
          "_gene_symbols_without_hgnc_id_mapping.tsv")))

    return(res)
  }
)


ds_hgnc_etbl_list <- imap(
  ds_hgnc_annot_list,
  function(hgnc_tbl, xname) {
    edge_tbl <- hgnc_tbl %>%
      select(cb_group_code, hgnc_code) %>%
      filter(!is.na(hgnc_code)) %>%
      mutate(predicate = "expressed_in") %>%
      rename(
        subject = cb_group_code,
        object = hgnc_code) %>%
      select(subject, predicate, object)

    stopifnot(identical(0L, sum(is.na(edge_tbl))))
    stopifnot(all(
      pull(edge_tbl, object) %in%
        pull(mapping_gsb_hgnc_tbl, hgnc_id)
    ))

    return(edge_tbl)
  }
)






out_node_tbl <- bind_rows(ds_node_tbl_list)
stopifnot(identical(0L, sum(is.na(out_node_tbl))))
stopifnot(identical(
  pull(out_node_tbl, node_id),
  unique(pull(out_node_tbl, node_id))))
stopifnot(identical(
  colnames(out_node_tbl),
  c("node_id", "node_label", "node_definition", "node_synonyms",
    "node_dbxrefs", "node_namespace")
))




ds_edge_tbl_list <- map(
  names(ds_node_tbl_list),
  function(xname) {
    edge_tbl_list <- map(
      list(
        ds_norm_uc_etbl_list,
        ds_umap_etbl_list,
        ds_hgnc_etbl_list
      ),
      function(etl) {
        return(etl[[xname]])
      }
    )


    walk(edge_tbl_list, function(e_tbl) {
      stopifnot(identical(0L, sum(is.na(e_tbl))))
      stopifnot(identical(
        unique(pull(e_tbl, subject)),
        pull(e_tbl, subject)
      ))
      stopifnot(identical(e_tbl, distinct(e_tbl)))
    })


    m_edge_tbl <- bind_rows(edge_tbl_list)
    stopifnot(identical(0L, sum(is.na(m_edge_tbl))))
    stopifnot(identical(m_edge_tbl, distinct(m_edge_tbl)))

    return(m_edge_tbl)
  }
)


out_edge_tbl <- bind_rows(ds_edge_tbl_list)
stopifnot(identical(0L, sum(is.na(out_edge_tbl))))
stopifnot(identical(
  distinct(out_edge_tbl),
  out_edge_tbl))
stopifnot(identical(
  colnames(out_edge_tbl),
  c("subject", "predicate", "object")
))
stopifnot(all(
  pull(out_edge_tbl, subject) %in%
    pull(out_node_tbl, node_id)))


out_sab_dir <- file.path(out_dir, SC_RNA_CLUST_SAB)
if (!dir.exists(out_sab_dir)) {
  dir.create(out_sab_dir)
}

write_tsv(
  out_node_tbl,
  file.path(out_sab_dir, "OWLNETS_node_metadata.txt"))

write_tsv(
  out_edge_tbl,
  file.path(out_sab_dir, "OWLNETS_edgelist.txt"))
