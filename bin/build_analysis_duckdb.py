#!/usr/bin/env python3
"""Build the QC-pass-only local analysis DuckDB."""

import argparse
import hashlib
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional, Union


EXPECTED_QC_PASS_GENOMES = 7_665


def parse_args() -> argparse.Namespace:
    repo_root = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=repo_root)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--metadata", type=Path, default=None)
    parser.add_argument("--pairwise-dir", type=Path, default=None)
    parser.add_argument("--markers", type=Path, default=None)
    parser.add_argument("--enterobacteriaceae", type=Path, default=None)
    parser.add_argument("--duckdb-bin", type=Path, default=None)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--smoke-rows", type=int, default=None)
    return parser.parse_args()


def sql_string(value: Union[str, Path]) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def find_duckdb(explicit: Optional[Path]) -> Path:
    candidates = []  # type: List[Path]
    if explicit is not None:
        candidates.append(explicit.expanduser())
    env_path = os.environ.get("DUCKDB_BIN")
    if env_path:
        candidates.append(Path(env_path).expanduser())
    path_match = shutil.which("duckdb")
    if path_match:
        candidates.append(Path(path_match))
    candidates.append(Path.home() / ".local/bin/duckdb")
    for candidate in candidates:
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate.resolve()
    raise FileNotFoundError("DuckDB CLI not found; supply --duckdb-bin or set DUCKDB_BIN")


def default_pairwise_dir(repo_root: Path) -> Path:
    local = repo_root / "data" / "all_tables_processed"
    if local.is_dir():
        return local
    return repo_root.parents[1] / repo_root.name / "data" / "all_tables_processed"


def sha256sum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def manifest_values(repo_root: Path, inputs: List[Path]) -> str:
    rows = []
    for path in inputs:
        try:
            display_path = path.relative_to(repo_root).as_posix()
        except ValueError:
            display_path = os.path.relpath(path, repo_root)
        rows.append(
            "(" + ", ".join(
                (sql_string(display_path), sql_string(sha256sum(path)), str(path.stat().st_size))
            ) + ")"
        )
    return ",\n        ".join(rows)


def build_sql(
    repo_root: Path,
    pair_glob: str,
    metadata_path: Path,
    marker_path: Path,
    enterobacteriaceae_path: Path,
    inputs: List[Path],
    threads: int,
    smoke_rows: Optional[int],
) -> str:
    row_limit = f"LIMIT {smoke_rows}" if smoke_rows is not None else ""
    mode = "smoke" if smoke_rows is not None else "production"
    built_at = datetime.now(timezone.utc).isoformat()
    manifest = manifest_values(repo_root, inputs)

    return f"""
SET threads = {threads};
SET preserve_insertion_order = false;

CREATE TABLE source_manifest (
    source_path VARCHAR NOT NULL,
    sha256 VARCHAR NOT NULL,
    size_bytes UBIGINT NOT NULL
);
INSERT INTO source_manifest VALUES
        {manifest};

CREATE TABLE build_metadata (
    key VARCHAR PRIMARY KEY,
    value VARCHAR NOT NULL
);
INSERT INTO build_metadata VALUES
    ('schema_version', '2'),
    ('build_mode', {sql_string(mode)}),
    ('built_at_utc', {sql_string(built_at)}),
    ('builder', 'bin/build_analysis_duckdb.py'),
    ('qc_source', 'assets/genome_tax_metadata.parquet'),
    ('qc_predicate', 'qc = ''pass''');

CREATE TEMP TABLE metadata_all AS
SELECT * FROM read_parquet({sql_string(metadata_path)});

SELECT CASE
    WHEN count(*) = count(DISTINCT ncbi_genome_accession) THEN count(*)
    ELSE error('genome_tax_metadata.parquet contains duplicate versioned accessions')
END AS validated_metadata_genomes
FROM metadata_all;

SELECT CASE
    WHEN count(*) = {EXPECTED_QC_PASS_GENOMES} THEN count(*)
    ELSE error('Unexpected number of QC-pass genomes')
END AS validated_qc_pass_genomes
FROM metadata_all
WHERE qc = 'pass';

CREATE TEMP VIEW metadata_accession_map AS
SELECT ncbi_genome_accession,
       split_part(ncbi_genome_accession, '.', 1) AS accession_base
FROM metadata_all;

SELECT CASE
    WHEN count(*) = count(DISTINCT accession_base) THEN count(*)
    ELSE error('Version-stripped accessions are not unique')
END AS validated_accession_bases
FROM metadata_accession_map;

CREATE TABLE genomes AS
SELECT row_number() OVER (ORDER BY ncbi_genome_accession)::INTEGER AS genome_id, *
FROM metadata_all
WHERE qc = 'pass'
ORDER BY ncbi_genome_accession;

CREATE TEMP VIEW accession_map AS
SELECT genome_id, ncbi_genome_accession,
       split_part(ncbi_genome_accession, '.', 1) AS accession_base
FROM genomes;

CREATE TEMP TABLE source_pairs AS
SELECT * FROM read_parquet({sql_string(pair_glob)}, union_by_name = true)
{row_limit};

SELECT CASE
    WHEN count(*) = 0 THEN count(*)
    ELSE error('At least one legacy pair accession is absent from canonical metadata')
END AS unmatched_pair_accessions
FROM (
    SELECT accession_base
    FROM (
        SELECT Genome1 AS accession_base FROM source_pairs
        UNION
        SELECT Genome2 AS accession_base FROM source_pairs
    ) source_accessions
    LEFT JOIN metadata_accession_map USING (accession_base)
    WHERE metadata_accession_map.accession_base IS NULL
);

CREATE TABLE pairwise_metrics AS
SELECT
    least(g1.genome_id, g2.genome_id)::INTEGER AS genome1_id,
    greatest(g1.genome_id, g2.genome_id)::INTEGER AS genome2_id,
    s.ANI::DOUBLE AS ani,
    s.AF::DOUBLE AS af,
    s.AAI::DOUBLE AS aai,
    s.SHARED_GENE_CONTENT::DOUBLE AS shared_gene_content,
    s.IDENT_18S::DOUBLE AS identity_18s,
    s.IDENT_28S::DOUBLE AS identity_28s,
    CASE
        WHEN m1.species_2026_09_01 IS NOT NULL
         AND m1.species_2026_09_01 = m2.species_2026_09_01 THEN 'species'
        WHEN m1.genus_2026_09_01 IS NOT NULL
         AND m1.genus_2026_09_01 = m2.genus_2026_09_01 THEN 'genus'
        WHEN m1.family_2026_09_01 IS NOT NULL
         AND m1.family_2026_09_01 = m2.family_2026_09_01 THEN 'family'
        WHEN m1.order_2026_09_01 IS NOT NULL
         AND m1.order_2026_09_01 = m2.order_2026_09_01 THEN 'order'
        WHEN m1.class_2026_09_01 IS NOT NULL
         AND m1.class_2026_09_01 = m2.class_2026_09_01 THEN 'class'
        WHEN m1.phylum_2026_09_01 IS NOT NULL
         AND m1.phylum_2026_09_01 = m2.phylum_2026_09_01 THEN 'phylum'
        WHEN m1.domain_2026_09_01 IS NOT NULL
         AND m1.domain_2026_09_01 = m2.domain_2026_09_01 THEN 'domain'
        ELSE NULL
    END AS lca_2026_09_01,
    lower(s.LSTR)::VARCHAR AS legacy_lstr
FROM source_pairs s
JOIN accession_map g1 ON s.Genome1 = g1.accession_base
JOIN accession_map g2 ON s.Genome2 = g2.accession_base
JOIN genomes m1 ON g1.genome_id = m1.genome_id
JOIN genomes m2 ON g2.genome_id = m2.genome_id
WHERE g1.genome_id <> g2.genome_id
ORDER BY genome1_id, genome2_id;

SELECT CASE
    WHEN count(*) = count(DISTINCT (genome1_id, genome2_id)) THEN count(*)
    ELSE error('Canonical pair identifiers are not unique')
END AS validated_pairs
FROM pairwise_metrics;

CREATE TEMP TABLE source_markers AS
SELECT * FROM read_parquet({sql_string(marker_path)})
{row_limit};

SELECT count(*) AS marker_accessions_not_in_metadata
FROM (
    SELECT accession
    FROM (
        SELECT genome1 AS accession FROM source_markers
        UNION
        SELECT genome2 AS accession FROM source_markers
    ) source_accessions
    LEFT JOIN metadata_all ON source_accessions.accession = metadata_all.ncbi_genome_accession
    WHERE metadata_all.ncbi_genome_accession IS NULL
);

CREATE TABLE marker_pairwise_results AS
SELECT
    least(g1.genome_id, g2.genome_id)::INTEGER AS genome1_id,
    greatest(g1.genome_id, g2.genome_id)::INTEGER AS genome2_id,
    m.*
FROM source_markers m
JOIN genomes g1 ON m.genome1 = g1.ncbi_genome_accession
JOIN genomes g2 ON m.genome2 = g2.ncbi_genome_accession
WHERE g1.genome_id <> g2.genome_id
ORDER BY genome1_id, genome2_id, region;

SELECT CASE
    WHEN count(*) = count(DISTINCT (genome1_id, genome2_id, region)) THEN count(*)
    ELSE error('Marker pair-region identifiers are not unique')
END AS validated_marker_pairs
FROM marker_pairwise_results;

CREATE TABLE enterobacteriaceae_fastani AS
SELECT * FROM read_parquet({sql_string(enterobacteriaceae_path)});

CREATE VIEW analysis_genomes AS SELECT * FROM genomes;
CREATE VIEW analysis_pairwise_metrics AS SELECT * FROM pairwise_metrics;

CREATE VIEW pairwise_accessions AS
SELECT g1.ncbi_genome_accession AS genome1,
       g2.ncbi_genome_accession AS genome2,
       p.* EXCLUDE (genome1_id, genome2_id)
FROM pairwise_metrics p
JOIN genomes g1 ON p.genome1_id = g1.genome_id
JOIN genomes g2 ON p.genome2_id = g2.genome_id;

SELECT CASE
    WHEN count(*) = 0 THEN count(*)
    ELSE error('Self comparisons remain in pairwise tables')
END AS validated_no_self_comparisons
FROM (
    SELECT genome1_id, genome2_id FROM pairwise_metrics WHERE genome1_id = genome2_id
    UNION ALL
    SELECT genome1_id, genome2_id FROM marker_pairwise_results WHERE genome1_id = genome2_id
);

INSERT INTO build_metadata
SELECT 'genome_rows', count(*)::VARCHAR FROM genomes
UNION ALL
SELECT 'pairwise_source_rows', count(*)::VARCHAR FROM source_pairs
UNION ALL
SELECT 'pairwise_rows', count(*)::VARCHAR FROM pairwise_metrics
UNION ALL
SELECT 'pairwise_excluded_rows',
       ((SELECT count(*) FROM source_pairs) - count(*))::VARCHAR
FROM pairwise_metrics
UNION ALL
SELECT 'marker_source_rows', count(*)::VARCHAR FROM source_markers
UNION ALL
SELECT 'marker_rows', count(*)::VARCHAR FROM marker_pairwise_results
UNION ALL
SELECT 'marker_excluded_rows',
       ((SELECT count(*) FROM source_markers) - count(*))::VARCHAR
FROM marker_pairwise_results
UNION ALL
SELECT 'enterobacteriaceae_rows', count(*)::VARCHAR FROM enterobacteriaceae_fastani;

ANALYZE;
CHECKPOINT;
"""


def main() -> int:
    args = parse_args()
    if args.threads < 1:
        raise ValueError("--threads must be at least 1")
    if args.smoke_rows is not None and args.smoke_rows < 1:
        raise ValueError("--smoke-rows must be at least 1")

    repo_root = args.repo_root.expanduser().resolve()
    pair_dir = args.pairwise_dir.expanduser().resolve() if args.pairwise_dir else default_pairwise_dir(repo_root).resolve()
    metadata_path = args.metadata.expanduser().resolve() if args.metadata else repo_root / "assets/genome_tax_metadata.parquet"
    marker_path = args.markers.expanduser().resolve() if args.markers else repo_root / "tmp_data/marker-pairwise-results.parquet"
    enterobacteriaceae_path = args.enterobacteriaceae.expanduser().resolve() if args.enterobacteriaceae else repo_root / "data/enterobacteriaceae_fastani.parquet"
    pair_files = sorted(pair_dir.glob("*.parquet"))
    inputs = [*pair_files, metadata_path, marker_path, enterobacteriaceae_path]
    missing = [path for path in inputs if not path.is_file()]
    if not pair_files or missing:
        details = "\n".join(str(path) for path in missing)
        raise FileNotFoundError(f"Missing required input files:\n{details}")

    output = (args.output or repo_root / "local_data/ani_microbial_eukaryotes.duckdb").expanduser().resolve()
    if output.exists() and not args.force:
        raise FileExistsError(f"Output already exists: {output}\nUse --force to replace it")
    output.parent.mkdir(parents=True, exist_ok=True)

    duckdb_bin = find_duckdb(args.duckdb_bin)
    temp_output = output.with_name(f".{output.name}.tmp-{os.getpid()}")
    temp_wal = Path(str(temp_output) + ".wal")
    for temporary in (temp_output, temp_wal):
        if temporary.exists():
            temporary.unlink()

    sql = build_sql(
        repo_root, (pair_dir / "*.parquet").as_posix(), metadata_path,
        marker_path, enterobacteriaceae_path, inputs, args.threads, args.smoke_rows,
    )
    print(f"Building {output}", flush=True)
    print(f"DuckDB CLI: {duckdb_bin}", flush=True)
    print(f"Metadata: {metadata_path}", flush=True)
    print(f"Pairwise data: {pair_dir}", flush=True)
    print(f"Markers: {marker_path}", flush=True)
    try:
        subprocess.run([str(duckdb_bin), str(temp_output), "-batch", "-c", sql], check=True)
        os.replace(temp_output, output)
    except BaseException:
        for temporary in (temp_output, temp_wal):
            if temporary.exists():
                temporary.unlink()
        raise

    size_mib = output.stat().st_size / (1024 * 1024)
    print(f"Created {output} ({size_mib:,.1f} MiB)")
    print(f"Inspect with:\n  {duckdb_bin} {output} -c \"SHOW TABLES;\"")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, FileExistsError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(2)
