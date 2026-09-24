#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
  library(igraph)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1]])), ".."))
metadata_path <- file.path(root, "data", "genome_tax_metadata.parquet")
matrix_dir <- file.path(root, "data")
if (!file.exists(file.path(matrix_dir, "fastani-upper-triangle.parquet"))) {
  matrix_dir <- file.path(root, "..", "..", basename(root), "data")
}
ani_path <- file.path(matrix_dir, "fastani-upper-triangle.parquet")
af_path <- file.path(matrix_dir, "fastani-AF-upper-triangle.parquet")
qpath <- function(path) gsub("'", "''", normalizePath(path, mustWork = TRUE), fixed = TRUE)

con <- dbConnect(duckdb(), dbdir = ":memory:")
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
invisible(dbExecute(con, "SET threads=4"))
invisible(dbExecute(con, "SET max_expression_depth=100000"))

vertices <- dbGetQuery(con, sprintf(
  "WITH base AS (
     SELECT ncbi_genome_accession AS name, species_2026_09_01 AS species
     FROM read_parquet('%s')
     WHERE qc = 'pass'
       AND species_2026_09_01 IS NOT NULL
       AND trim(species_2026_09_01) <> ''
   ), eligible AS (
     SELECT species FROM base GROUP BY species HAVING count(*) >= 3
   )
   SELECT b.name, b.species FROM base b JOIN eligible e USING (species) ORDER BY b.name",
  qpath(metadata_path)
))

if (nrow(vertices) == 0) stop("No eligible genomes", call. = FALSE)

duckdb_register(con, "eligible_vertices", vertices)
edges <- dbGetQuery(con, sprintf(
  "WITH ani_long AS (
     SELECT assembly_accession AS source, compared_accession AS target, value::DOUBLE AS ani
     FROM (
       UNPIVOT read_parquet('%s')
       ON COLUMNS(* EXCLUDE (assembly_accession))
       INTO NAME compared_accession VALUE value
     )
     WHERE assembly_accession <> compared_accession
       AND value IS NOT NULL AND isfinite(value::DOUBLE) AND value::DOUBLE >= 95
   ), af_long AS (
     SELECT assembly_accession AS source, compared_accession AS target, value::DOUBLE AS af
     FROM (
       UNPIVOT read_parquet('%s')
       ON COLUMNS(* EXCLUDE (assembly_accession))
       INTO NAME compared_accession VALUE value
     )
     WHERE assembly_accession <> compared_accession
       AND value IS NOT NULL AND isfinite(value::DOUBLE) AND value::DOUBLE >= 0.6
   )
   SELECT DISTINCT least(a.source, a.target) AS source,
                   greatest(a.source, a.target) AS target
   FROM ani_long a
   JOIN af_long f USING (source, target)
   JOIN eligible_vertices v1 ON a.source = v1.name
   JOIN eligible_vertices v2 ON a.target = v2.name",
  qpath(ani_path), qpath(af_path)
))
duckdb_unregister(con, "eligible_vertices")

g <- graph_from_data_frame(edges, directed = FALSE, vertices = vertices)
g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
cc <- components(g)
membership <- data.frame(
  name = V(g)$name,
  species = V(g)$species,
  component = as.integer(cc$membership),
  stringsAsFactors = FALSE
)

component_species <- unique(membership[c("component", "species")])
n_species <- table(component_species$component)
component_ids <- seq_len(cc$no)
component_species_counts <- as.integer(n_species[as.character(component_ids)])
component_species_counts[is.na(component_species_counts)] <- 0L
homogeneous_ids <- component_ids[component_species_counts == 1L]
heterogeneous_ids <- component_ids[component_species_counts > 1L]

species_components <- split(component_species$component, component_species$species)
recovery <- vapply(species_components, function(ids) {
  ids <- unique(ids)
  if (any(component_species_counts[ids] > 1L)) {
    "heterogeneous"
  } else if (length(ids) == 1L) {
    "single_homogeneous"
  } else {
    "multiple_homogeneous"
  }
}, character(1))

outliers <- vapply(heterogeneous_ids, function(id) {
  counts <- table(membership$species[membership$component == id])
  as.integer(sum(counts) - max(counts))
}, integer(1))

clique_sizes <- integer()
clique_species_counts <- integer()
for (component_id in component_ids) {
  sg <- induced_subgraph(g, which(cc$membership == component_id))
  cliques <- max_cliques(sg, min = 2)
  if (length(cliques) == 0) next
  clique_sizes <- c(clique_sizes, vapply(cliques, length, integer(1)))
  clique_species_counts <- c(
    clique_species_counts,
    vapply(cliques, function(clique) {
      length(unique(vertex_attr(sg, "species", index = clique)))
    }, integer(1))
  )
}

stat <- function(metric, value, denominator = NA_real_) {
  data.frame(
    metric = metric,
    value = as.numeric(value),
    denominator = as.numeric(denominator),
    percent = if (!is.na(denominator) && denominator > 0) 100 * value / denominator else NA_real_
  )
}

n_species_total <- length(species_components)
n_components <- cc$no
n_heterogeneous <- length(heterogeneous_ids)
n_cliques <- length(clique_sizes)
stats <- do.call(rbind, list(
  stat("eligible QC-pass genomes", vcount(g)),
  stat("eligible resolved species with at least 3 representatives", n_species_total),
  stat("qualifying ANI>=95 and AF>=0.6 edges", ecount(g)),
  stat("connected components", n_components),
  stat("isolated genomes", sum(sizes(cc) == 1L), vcount(g)),
  stat("homogeneous connected components", length(homogeneous_ids), n_components),
  stat("heterogeneous connected components", n_heterogeneous, n_components),
  stat("species in one homogeneous component", sum(recovery == "single_homogeneous"), n_species_total),
  stat("species in multiple homogeneous components", sum(recovery == "multiple_homogeneous"), n_species_total),
  stat("species in at least one heterogeneous component", sum(recovery == "heterogeneous"), n_species_total),
  stat("heterogeneous components with exactly 1 outlier genome", sum(outliers == 1L), n_heterogeneous),
  stat("heterogeneous components with <=5 outlier genomes", sum(outliers <= 5L), n_heterogeneous),
  stat("maximal cliques with at least 2 genomes", n_cliques),
  stat("homogeneous maximal cliques with at least 2 genomes", sum(clique_species_counts == 1L), n_cliques),
  stat("heterogeneous maximal cliques with at least 2 genomes", sum(clique_species_counts > 1L), n_cliques)
))

write.csv(stats, file = "", row.names = FALSE, na = "")
