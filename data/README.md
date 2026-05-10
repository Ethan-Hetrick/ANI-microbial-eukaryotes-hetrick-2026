# Data

`all_tables_processed/` is a five-part Parquet dataset containing the processed
pairwise comparison table used by the publication notebook. The table is split
only to keep each individual file below GitHub's 100 MB file-size limit; the
parts do not need to be concatenated before use.

Summary:

- Rows: 25,156,357
- Columns: 27
- Format: Parquet, ZSTD compression

Column reference:

- Genome pair identifiers: `Genome1`, `Genome2`
- Whole-genome nucleotide similarity: `ANI`, average nucleotide identity; `AF`, FastANI alignment fraction.
- Amino-acid/gene-content similarity: `AAI`, average amino acid identity; `SHARED_GENE_CONTENT`, fraction of shared protein-coding genes.
- rRNA identity: `IDENT_18S`, 18S rRNA nucleotide identity; `IDENT_28S`, 28S rRNA nucleotide identity.
- Genome 1 taxonomy: `GEN1_DOMAIN`, `GEN1_PHYLUM`, `GEN1_CLASS`, `GEN1_ORDER`, `GEN1_FAMILY`, `GEN1_GENUS`, `GEN1_SPECIES`.
- Genome 2 taxonomy: `GEN2_DOMAIN`, `GEN2_PHYLUM`, `GEN2_CLASS`, `GEN2_ORDER`, `GEN2_FAMILY`, `GEN2_GENUS`, `GEN2_SPECIES`.
- Assembly quality summaries: `GEN1_N50`, `GEN1_L50`, `GEN2_N50`, `GEN2_L50`.
- Pairwise taxonomic relationship: `LSTR`, lowest shared taxonomic rank for the genome pair.

Additional notebook support file in this directory:

- `enterobacteriaceae_fastani.parquet`: FastANI data for the eukaryote-vs-prokaryote comparison.

`enterobacteriaceae_fastani.parquet` columns:

- `query_accession`, `reference_accession`: genome accessions in each FastANI comparison.
- `ANI`: average nucleotide identity.
- `mapped_frags`, `total_frags`: mapped and total FastANI fragments.
- `AF`: alignment fraction, computed as `mapped_frags / total_frags`.

The Figure 1 genome taxonomy and genus annotation tables are stored under
`../assets/`.

Read the split dataset directly with PyArrow in Python or DuckDB in R, for
example `data/all_tables_processed/*.parquet`.
