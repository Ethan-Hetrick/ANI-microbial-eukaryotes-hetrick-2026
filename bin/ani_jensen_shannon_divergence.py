#!/usr/bin/env python3

from pathlib import Path

import duckdb
import numpy as np
from scipy.stats import entropy


repo_root = Path(__file__).resolve().parents[1]
database = repo_root / "local_data" / "ani_microbial_eukaryotes.duckdb"

with duckdb.connect(str(database), read_only=True) as con:
    data = con.execute(
        """
        SELECT ani, lca_2026_09_01 = 'species' AS intraspecies
        FROM pairwise_metrics
        WHERE ani IS NOT NULL
          AND isfinite(ani)
          AND lca_2026_09_01 IS NOT NULL
        """
    ).fetchnumpy()

ani = np.asarray(data["ani"])
intraspecies = np.asarray(data["intraspecies"])
bin_width = 0.1
lower_bound = np.floor(ani.min() / bin_width) * bin_width
upper_bound = np.ceil(ani.max() / bin_width) * bin_width
bin_count = round((upper_bound - lower_bound) / bin_width)
bins = lower_bound + np.arange(bin_count + 1) * bin_width
bins[-1] = np.nextafter(upper_bound, np.inf)

intra_counts = np.histogram(ani[intraspecies], bins=bins)[0].astype(float)
inter_counts = np.histogram(ani[~intraspecies], bins=bins)[0].astype(float)

intra_probability = (intra_counts + 0.5) / (intra_counts.sum() + 0.5 * len(intra_counts))
inter_probability = (inter_counts + 0.5) / (inter_counts.sum() + 0.5 * len(inter_counts))
mixture = 0.5 * (intra_probability + inter_probability)
js_divergence = 0.5 * entropy(intra_probability, mixture) + 0.5 * entropy(inter_probability, mixture)

print(f"intraspecies_comparisons\t{int(intra_counts.sum())}")
print(f"interspecies_comparisons\t{int(inter_counts.sum())}")
print(f"bin_width_percent_identity\t{bin_width}")
print(f"Jensen_Shannon_divergence_nats\t{js_divergence:.6f}")
