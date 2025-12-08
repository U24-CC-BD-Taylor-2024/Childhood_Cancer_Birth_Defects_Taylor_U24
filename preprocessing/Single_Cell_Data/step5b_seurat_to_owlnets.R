library(tidyverse)



analysis_id <- paste0(
  "seurat_to_owlnets_v1"
)
out_dir <- paste0("results/", analysis_id, "_results")
out_int_dir <- paste0("results/", analysis_id, "_intermediate_files")






node_files <- dir(
  out_dir, pattern = "OWLNETS_node_metadata\\.txt",
  recursive = TRUE)

extract_sab <- function(code) {
  code_match <- str_match(code, "^([^:]+):[^ :_]+$")
  stopifnot(identical(0L, sum(is.na(code_match))))

  res <- code_match[, 2]
  stopifnot(is.character(res))

  return(res)
}

node_stats_tbl <- bind_rows(map(node_files, function(x) {
  xtbl <- read_tsv(
    file.path(out_dir, x), na = "NA", col_types = "cccccc",
    progress = FALSE)

  stopifnot(identical(0L, sum(is.na(xtbl))))

  stopifnot(identical(
    pull(xtbl, node_id), unique(pull(xtbl, node_id))))

  sab_tbl <- xtbl %>%
    mutate(SAB = extract_sab(node_id))

  stopifnot(identical(0L, sum(is.na(sab_tbl))))

  stats_tbl <- sab_tbl %>%
    count(SAB, name = "n_nodes") %>%
    mutate(file = x) %>%
    relocate(file) %>%
    arrange(file, SAB)

  stopifnot(identical(0L, sum(is.na(stats_tbl))))

  return(stats_tbl)
}))



edge_files <- dir(
  out_dir, pattern = "OWLNETS_edgelist\\.txt",
  recursive = TRUE)

edge_stats_tbl <- bind_rows(map(edge_files, function(x) {
  xtbl <- read_tsv(
    file.path(out_dir, x), na = "NA", col_types = "ccc",
    progress = FALSE)

  stopifnot(identical(0L, sum(is.na(xtbl))))
  
  all_sab_stats_tbl <- tibble(
    file = x,
    subject_SAB = "All_subject_SAB",
    predicate = "All_predicate",
    object_SAB = "All_object_SAB",
    n_edges = nrow(xtbl)
  )

  if (nrow(xtbl) == 0) {
    res <- all_sab_stats_tbl
  } else {
    stopifnot(nrow(xtbl) > 0)
    stopifnot(identical(xtbl[, 1:ncol(xtbl)], distinct(xtbl)))

    sab_tbl <- xtbl %>%
      mutate(
        subject_SAB = extract_sab(subject),
        object_SAB = extract_sab(object)
      )

    stopifnot(identical(0L, sum(is.na(sab_tbl))))

    stats_tbl <- sab_tbl %>%
      count(
        subject_SAB, predicate, object_SAB,
        name = "n_edges") %>%
      mutate(file = x) %>%
      relocate(file) %>%
      arrange(file, subject_SAB, predicate, object_SAB)

    stopifnot(identical(0L, sum(is.na(stats_tbl))))
    stopifnot(identical(
      pull(all_sab_stats_tbl, n_edges),
      sum(pull(stats_tbl, n_edges))))

    res <- bind_rows(all_sab_stats_tbl, stats_tbl)
  }

  stopifnot(identical(0L, sum(is.na(res))))

  return(res)
}))



write_tsv(
  node_stats_tbl,
  file.path(
    out_dir,
    paste0("scRNA_seq", "_node_stats.tsv")))

write_tsv(
  edge_stats_tbl,
  file.path(
    out_dir,
    paste0("scRNA_seq", "_edge_stats.tsv")))




save_files <- dir(out_dir, recursive = TRUE)
save_files <- save_files[!(save_files == "sha256sum.txt")]
print(save_files)

save_file_sha256 <- tools::sha256sum(
  file.path(out_dir, save_files))

save_file_sha256_line <- paste(
  save_file_sha256, names(save_file_sha256),
  sep = "  ")

write_lines(
  save_file_sha256_line,
  file.path(out_dir, "sha256sum.txt"))
