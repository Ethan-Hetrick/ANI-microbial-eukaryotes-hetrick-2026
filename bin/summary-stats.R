#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1]])), ".."))
database_path <- normalizePath(file.path(root, "local_data", "ani_microbial_eukaryotes.duckdb"))

con <- dbConnect(duckdb(), dbdir = database_path, read_only = TRUE)
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

counts <- dbGetQuery(con,
  "WITH pairs AS (
     SELECT p.ani,
            g1.species_2026_09_01 = g2.species_2026_09_01 AS same_species
     FROM pairwise_metrics p
     JOIN genomes g1 ON p.genome1_id = g1.genome_id
     JOIN genomes g2 ON p.genome2_id = g2.genome_id
     WHERE p.ani IS NOT NULL
       AND isfinite(p.ani)
       AND g1.species_2026_09_01 IS NOT NULL
       AND trim(g1.species_2026_09_01) <> ''
       AND g2.species_2026_09_01 IS NOT NULL
       AND trim(g2.species_2026_09_01) <> ''
   )
   SELECT
     count(*) AS total_comparisons,
     count(*) FILTER (WHERE same_species) AS same_species_comparisons,
     count(*) FILTER (WHERE NOT same_species) AS different_species_comparisons,
     count(*) FILTER (WHERE same_species AND ani >= 95) AS true_positive,
     count(*) FILTER (WHERE same_species AND ani < 95) AS false_negative,
     count(*) FILTER (WHERE NOT same_species AND ani >= 95) AS false_positive,
     count(*) FILTER (WHERE NOT same_species AND ani < 95) AS true_negative
   FROM pairs")

for (name in names(counts)) counts[[name]] <- as.numeric(counts[[name]])
counts$same_species_percent <- 100 * counts$same_species_comparisons / counts$total_comparisons
counts$different_species_percent <- 100 * counts$different_species_comparisons / counts$total_comparisons
counts$recall_percent <- 100 * counts$true_positive / counts$same_species_comparisons
counts$false_negative_percent <- 100 * counts$false_negative / counts$same_species_comparisons
counts$false_positive_rate_percent <- 100 * counts$false_positive / counts$different_species_comparisons
counts$precision_percent <- 100 * counts$true_positive / (counts$true_positive + counts$false_positive)
counts$specificity_percent <- 100 * counts$true_negative / (counts$true_negative + counts$false_positive)

write.csv(counts, file = "", row.names = FALSE, na = "")
