library(DBI)
library(duckdb)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1]])), ".."))
database_path <- normalizePath(file.path(root, "local_data", "ani_microbial_eukaryotes.duckdb"))
metadata_path <- normalizePath(file.path(root, "assets", "genome_tax_metadata.parquet"))
qpath <- function(path) gsub("'", "''", path, fixed = TRUE)

con <- dbConnect(
  duckdb(),
  dbdir = database_path,
  read_only = TRUE
)

pairs <- dbGetQuery(con, sprintf("
  WITH retained AS (
    SELECT ncbi_genome_accession AS accession,
           genus_2026_09_01 AS genus,
           species_2026_09_01 AS species
    FROM read_parquet('%s')
    WHERE qc = 'pass'
  )
  SELECT
    CASE
      WHEN m1.species = m2.species
        THEN 'same species'
      WHEN m1.genus = m2.genus
       AND m1.species <> m2.species
        THEN 'different species within genus'
    END AS comparison_type,
    p.ani
  FROM analysis_pairwise_metrics p
  JOIN genomes g1 ON p.genome1_id = g1.genome_id
  JOIN genomes g2 ON p.genome2_id = g2.genome_id
  JOIN retained m1 ON g1.ncbi_genome_accession = m1.accession
  JOIN retained m2 ON g2.ncbi_genome_accession = m2.accession
  WHERE p.genome1_id <> p.genome2_id
    AND p.ani IS NOT NULL
    AND isfinite(p.ani)
    AND (
      m1.species = m2.species
      OR (
        m1.genus = m2.genus
        AND m1.species <> m2.species
      )
    )
", qpath(metadata_path)))

dbDisconnect(con, shutdown = TRUE)

bootstrap_mean <- function(x, replicates = 1000, seed = 2026) {
  set.seed(seed)

  boot_means <- replicate(
    replicates,
    mean(x[sample.int(length(x), replace = TRUE)])
  )

  c(
    comparisons = length(x),
    mean = mean(x),
    ci_low = unname(quantile(boot_means, 0.025)),
    ci_high = unname(quantile(boot_means, 0.975))
  )
}

results <- do.call(
  rbind,
  lapply(
    split(pairs$ani, pairs$comparison_type),
    bootstrap_mean
  )
)

print(results)
