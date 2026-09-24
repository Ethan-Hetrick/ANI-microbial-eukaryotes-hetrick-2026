library(DBI)
library(duckdb)

con <- dbConnect(
  duckdb(),
  dbdir = "local_data/ani_microbial_eukaryotes.duckdb",
  read_only = TRUE
)

pairs <- dbGetQuery(con, "
  SELECT
    CASE
      WHEN g1.species_2026_09_01 = g2.species_2026_09_01
        THEN 'same species'
      WHEN g1.genus_2026_09_01 = g2.genus_2026_09_01
       AND g1.species_2026_09_01 <> g2.species_2026_09_01
        THEN 'different species within genus'
    END AS comparison_type,
    p.ani
  FROM analysis_pairwise_metrics p
  JOIN genomes g1 ON p.genome1_id = g1.genome_id
  JOIN genomes g2 ON p.genome2_id = g2.genome_id
  WHERE g1.qc = 'pass'
    AND g2.qc = 'pass'
    AND p.genome1_id <> p.genome2_id
    AND p.ani IS NOT NULL
    AND isfinite(p.ani)
    AND (
      g1.species_2026_09_01 = g2.species_2026_09_01
      OR (
        g1.genus_2026_09_01 = g2.genus_2026_09_01
        AND g1.species_2026_09_01 <> g2.species_2026_09_01
      )
    )
")

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

