#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
  library(igraph)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1]])), ".."))
database_path <- normalizePath(file.path(root, "local_data", "ani_microbial_eukaryotes.duckdb"))

con <- dbConnect(duckdb(), dbdir = database_path, read_only = TRUE)
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
invisible(dbExecute(con, "SET threads=4"))
invisible(dbExecute(con, "SET max_expression_depth=100000"))

vertices <- dbGetQuery(con,
  "WITH base AS (
     SELECT genome_id, ncbi_genome_accession AS name, species_2026_09_01 AS species
     FROM genomes
     WHERE species_2026_09_01 IS NOT NULL
       AND trim(species_2026_09_01) <> ''
   ), eligible AS (
     SELECT species FROM base GROUP BY species HAVING count(*) >= 3
   )
   SELECT b.genome_id, b.name, b.species
   FROM base b JOIN eligible e USING (species)
   ORDER BY b.name")

if (nrow(vertices) == 0) stop("No eligible genomes", call. = FALSE)

duckdb_register(con, "eligible_vertices", vertices)
edges <- dbGetQuery(con,
  "SELECT g1.ncbi_genome_accession AS source,
          g2.ncbi_genome_accession AS target
   FROM pairwise_metrics p
   JOIN eligible_vertices v1 ON p.genome1_id = v1.genome_id
   JOIN eligible_vertices v2 ON p.genome2_id = v2.genome_id
   JOIN genomes g1 ON p.genome1_id = g1.genome_id
   JOIN genomes g2 ON p.genome2_id = g2.genome_id
   WHERE p.ani >= 95 AND p.af >= 0.6")
duckdb_unregister(con, "eligible_vertices")

vertices$genome_id <- NULL

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
