#!/usr/bin/env Rscript
# Rebuild WHO annotation tables from data/who_compendium/current.xlsx and
# report any AA markers that lack genomic coordinates in the principal table.
#
# Usage (from repo root):
#   Rscript scripts/update_who_compendium.R
#
# Optional args:
#   --xlsx PATH     Use a specific Excel file instead of current.xlsx
#   --set-current   Also copy --xlsx to data/who_compendium/current.xlsx
#                   and archive it (requires --version LABEL)
#   --version LABEL Version label for archive / VERSION.txt (e.g. 2026-12_v1.1)

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[[i + 1]]
}

has_flag <- function(flag) flag %in% args

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop(
    "Package 'readxl' is required.\n",
    "Install it with: install.packages(\"readxl\")",
    call. = FALSE
  )
}
if (!requireNamespace("readr", quietly = TRUE) ||
    !requireNamespace("dplyr", quietly = TRUE) ||
    !requireNamespace("tidyr", quietly = TRUE) ||
    !requireNamespace("stringr", quietly = TRUE)) {
  stop(
    "Packages readr, dplyr, tidyr, and stringr are required.\n",
    "Install with: install.packages(c(\"readr\", \"dplyr\", \"tidyr\", \"stringr\"))",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
who_dir <- file.path(root, "data", "who_compendium")
archive_dir <- file.path(who_dir, "archive")
current_xlsx <- file.path(who_dir, "current.xlsx")
version_file <- file.path(who_dir, "VERSION.txt")
principal_path <- file.path(root, "data", "principal_resistance_marker_info_table.tsv")
long_out <- file.path(root, "data", "who_compendium_annotations_long.tsv")
summary_out <- file.path(root, "data", "who_compendium_marker_summary.tsv")
missing_out <- file.path(root, "data", "who_compendium_missing_coords.tsv")

xlsx_arg <- get_arg("--xlsx")
version_label <- get_arg("--version")
set_current <- has_flag("--set-current")

if (set_current) {
  if (is.null(xlsx_arg) || is.null(version_label)) {
    stop("--set-current requires both --xlsx PATH and --version LABEL", call. = FALSE)
  }
  if (!file.exists(xlsx_arg)) stop("Excel file not found: ", xlsx_arg, call. = FALSE)
  dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)
  archive_path <- file.path(archive_dir, paste0(version_label, ".xlsx"))
  file.copy(xlsx_arg, archive_path, overwrite = TRUE)
  file.copy(xlsx_arg, current_xlsx, overwrite = TRUE)
  writeLines(
    c(
      version_label,
      "Source: WHO Compendium of molecular markers for antimalarial drug resistance",
      paste("Archived as:", basename(archive_path)),
      "File retained in original Excel format."
    ),
    version_file
  )
  message("Set current WHO compendium to ", version_label)
  message("  current: ", current_xlsx)
  message("  archive: ", archive_path)
}

xlsx_path <- if (!is.null(xlsx_arg) && !set_current) {
  normalizePath(xlsx_arg, winslash = "/", mustWork = TRUE)
} else {
  if (!file.exists(current_xlsx)) {
    stop(
      "Missing ", current_xlsx, "\n",
      "Place the WHO Excel file there, or run with:\n",
      "  Rscript scripts/update_who_compendium.R --set-current --xlsx PATH --version LABEL",
      call. = FALSE
    )
  }
  current_xlsx
}

message("Reading WHO Excel: ", xlsx_path)
if (file.exists(version_file)) {
  message("VERSION.txt:\n  ", paste(readLines(version_file, warn = FALSE), collapse = "\n  "))
}

who_to_princ <- c(
  pfk13 = "k13", k13 = "k13", kelch13 = "k13",
  pfcrt = "crt", crt = "crt",
  pfmdr1 = "mdr1", mdr1 = "mdr1",
  pfdhfr = "dhfr", dhfr = "dhfr",
  pfdhps = "dhps", dhps = "dhps",
  pfcytb = "cytb", cytb = "cytb",
  pfcoronin = "coronin", coronin = "coronin",
  pfgch1 = "gch1", gch1 = "gch1",
  `pfpm2/3` = "plasmepsin2_3"
)

normalize_gene <- function(g) {
  key <- gsub("\\s+", "", tolower(trimws(as.character(g))))
  if (key %in% names(who_to_princ)) who_to_princ[[key]] else key
}

class_bucket <- function(c) {
  c <- tolower(as.character(c))
  dplyr::case_when(
    grepl("validated", c) ~ "Validated",
    grepl("candidate", c) ~ "Candidate",
    grepl("potential", c) ~ "Potential",
    TRUE ~ "Other"
  )
}

parse_aas <- function(text) {
  matches <- str_match_all(as.character(text), "([A-Za-z]?)([0-9]+)([A-Za-z](/[A-Za-z]*)?)?")[[1]]
  if (nrow(matches) == 0) {
    return(tibble(wt = character(), aa_position = integer(), mut = character()))
  }
  tibble(
    wt = ifelse(matches[, 2] == "", NA_character_, toupper(matches[, 2])),
    aa_position = as.integer(matches[, 3]),
    mut = ifelse(is.na(matches[, 4]) | matches[, 4] == "", NA_character_, toupper(matches[, 4]))
  ) %>%
    filter(!(is.na(.data$wt) & is.na(.data$mut)))
}

expand_row <- function(gene_who, alteration) {
  alteration <- as.character(alteration)
  gene_who <- as.character(gene_who)
  if (grepl("amplif", alteration, ignore.case = TRUE)) {
    return(tibble(
      gene = normalize_gene(gene_who),
      aa_position = NA_integer_,
      kind = "amplification"
    ))
  }

  # Dual-gene haplotype rows: "Pfdhfr & Pfdhps" with alteration split on &
  if (grepl("&", gene_who) && grepl("pfdhfr|pfdhps", gene_who, ignore.case = TRUE)) {
    parts <- trimws(strsplit(alteration, "&", fixed = TRUE)[[1]])
    genes <- vapply(trimws(strsplit(gene_who, "&", fixed = TRUE)[[1]]), normalize_gene, character(1))
    if (length(parts) == length(genes)) {
      out <- bind_rows(lapply(seq_along(parts), function(i) {
        aas <- parse_aas(parts[[i]])
        if (nrow(aas) == 0) {
          return(tibble(gene = genes[[i]], aa_position = NA_integer_, kind = "unparsed"))
        }
        mutate(aas, gene = genes[[i]], kind = "aa") %>%
          select(gene, aa_position, kind)
      }))
      return(out)
    }
    aas <- parse_aas(alteration)
    return(
      aas %>%
        mutate(
          gene = ifelse(.data$aa_position < 300, "dhfr", "dhps"),
          kind = "aa"
        ) %>%
        select(gene, aa_position, kind)
    )
  }

  gene <- normalize_gene(gene_who)
  aas <- parse_aas(alteration)
  if (nrow(aas) == 0) {
    return(tibble(gene = gene, aa_position = NA_integer_, kind = "unparsed"))
  }
  mutate(aas, gene = gene, kind = "aa") %>%
    select(gene, aa_position, kind)
}

skip_sheets <- c(
  "User guide",
  "Combined  (Pf)",
  "References  (Pf)",
  "PfCRT haplotypes (for ref.)"
)

sheets <- setdiff(excel_sheets(xlsx_path), skip_sheets)
frames <- lapply(sheets, function(sheet) {
  df <- read_excel(xlsx_path, sheet = sheet)
  names(df) <- str_replace_all(names(df), "\\n", " ")
  names(df) <- str_trim(names(df))
  if (!"Drug" %in% names(df)) {
    df$Drug <- str_trim(str_replace(sheet, " \\(Pf\\)$", ""))
  }
  needed <- c("Drug", "Gene(s)", "Alteration(s)", "Classification")
  missing_cols <- setdiff(needed, names(df))
  if (length(missing_cols) > 0) {
    warning("Skipping sheet ", sheet, "; missing columns: ", paste(missing_cols, collapse = ", "))
    return(NULL)
  }
  df %>%
    select(
      drug = `Drug`,
      gene_who = `Gene(s)`,
      alteration = `Alteration(s)`,
      classification = `Classification`
    ) %>%
    mutate(across(everything(), ~ str_trim(as.character(.x)))) %>%
    filter(
      !is.na(.data$gene_who), .data$gene_who != "", .data$gene_who != "Gene(s)",
      !is.na(.data$alteration), .data$alteration != "", .data$alteration != "nan"
    )
})

comb <- bind_rows(frames) %>%
  mutate(class_bucket = class_bucket(.data$classification))

expanded <- bind_rows(lapply(seq_len(nrow(comb)), function(i) {
  r <- comb[i, , drop = FALSE]
  expand_row(r$gene_who, r$alteration) %>%
    mutate(
      drug = r$drug,
      classification = r$classification,
      class_bucket = r$class_bucket,
      alteration = r$alteration,
      gene_who = r$gene_who
    )
}))

amps <- expanded %>%
  filter(.data$kind == "amplification") %>%
  distinct(.data$drug, .data$gene, .data$gene_who, .data$classification, .data$class_bucket)

long <- expanded %>%
  filter(.data$kind == "aa", !is.na(.data$aa_position)) %>%
  distinct(
    .data$gene, .data$aa_position, .data$drug, .data$classification,
    .data$class_bucket, .data$alteration, .data$gene_who
  ) %>%
  arrange(.data$gene, .data$aa_position, .data$drug, .data$class_bucket)

join_sorted <- function(x) paste(sort(unique(as.character(x))), collapse = "; ")

summary_tbl <- long %>%
  group_by(.data$gene, .data$aa_position) %>%
  summarise(
    who_drugs = join_sorted(.data$drug),
    who_class_buckets = join_sorted(.data$class_bucket),
    who_classifications = join_sorted(.data$classification),
    who_alterations = join_sorted(.data$alteration),
    .groups = "drop"
  ) %>%
  mutate(in_who_compendium = TRUE)

if (!file.exists(principal_path)) {
  stop("Missing principal marker table: ", principal_path, call. = FALSE)
}

principal <- read_tsv(principal_path, show_col_types = FALSE) %>%
  mutate(
    gene = as.character(.data$gene),
    aa_position = as.integer(.data$aa_position)
  )

princ_keys <- principal %>%
  distinct(.data$gene, .data$aa_position) %>%
  mutate(has_coords = TRUE)

missing_aa <- summary_tbl %>%
  mutate(aa_position = as.integer(.data$aa_position)) %>%
  left_join(princ_keys, by = c("gene", "aa_position")) %>%
  filter(is.na(.data$has_coords)) %>%
  select(-has_coords, -in_who_compendium) %>%
  arrange(.data$gene, .data$aa_position)

write_tsv(long, long_out)
write_tsv(summary_tbl, summary_out)
write_tsv(missing_aa, missing_out)

message("")
message("Wrote:")
message("  ", long_out, "  (", nrow(long), " rows)")
message("  ", summary_out, "  (", nrow(summary_tbl), " AA markers)")
message("  ", missing_out, "  (", nrow(missing_aa), " AA markers without coords)")

message("")
message("Coverage vs principal table:")
message("  WHO unique AA markers: ", nrow(summary_tbl))
message("  With genomic coords:   ", nrow(summary_tbl) - nrow(missing_aa))
message("  Missing coords:        ", nrow(missing_aa))

if (nrow(missing_aa) > 0) {
  message("")
  message("AA markers still needing genomic locations:")
  for (i in seq_len(nrow(missing_aa))) {
    r <- missing_aa[i, ]
    message(sprintf(
      "  %-12s aa %-4s  drugs=%s  class=%s",
      r$gene, r$aa_position, r$who_drugs, r$who_class_buckets
    ))
  }
  message("Details: ", missing_out)
} else {
  message("")
  message("All WHO AA markers have genomic coordinates in the principal table.")
}

if (nrow(amps) > 0) {
  message("")
  message("Amplification markers (not AA sites; optional separate handling):")
  for (i in seq_len(nrow(amps))) {
    r <- amps[i, ]
    message(sprintf(
      "  %-14s  %-22s  %s",
      r$gene, r$drug, r$class_bucket
    ))
  }
}

message("")
message("Done. Restart the Shiny app to load the refreshed WHO annotations.")
