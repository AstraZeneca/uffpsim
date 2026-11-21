"""Summarise per-query outcomes into a Markdown table (results/summary.md)."""
import csv
import sys
from collections import Counter, defaultdict

from common import EXPERIMENT_DIR, database_dir, load_config


def read_outcomes(path):
    counts, with_hits, max_error = defaultdict(Counter), defaultdict(set), 0.0
    with open(path) as f:
        for row in csv.DictReader(f):
            counts[(row["threshold"], row["k"], row["method"])][row["outcome"]] += 1
            if int(row["n_expected"]) > 0:
                with_hits[(row["threshold"], row["k"])].add(row["query"])
            max_error = max(max_error, float(row["max_score_error"]))
    return counts, with_hits, max_error


def cell(outcomes):
    total = sum(outcomes.values())
    text = f"{outcomes['identical']:,} / {total:,} ({100 * outcomes['identical'] / total:.2f}%)"
    deviations = [f"{n:,} {name}" for name, n in sorted(outcomes.items()) if name != "identical"]
    return text + (f" [{', '.join(deviations)}]" if deviations else "")


def table(config, database):
    counts, with_hits, max_error = read_outcomes(database_dir(config, database) / "per_query.csv")
    methods = [f"UFFPSim {mode}" for mode in config["uffpsim"]["modes"]] + ["FPSim2"]
    lines = [
        f"### {database['fp_size']}-bit fingerprints (inner-clustering threshold {database['inner_clustering_threshold']})",
        "",
        "| Search | Threshold | Queries with hits | " + " | ".join(methods) + " |",
        "|---|---|---|" + "---|" * len(methods),
    ]
    for threshold in config["searches"]["thresholds"]:
        for k in ["all", *config["searches"]["top_k"]]:
            search = "all hits" if k == "all" else f"top-{k}"
            cells = [cell(counts[(str(threshold), str(k), method)]) for method in methods]
            n_with_hits = len(with_hits[(str(threshold), str(k))])
            lines.append(f"| {search} | {threshold} | {n_with_hits:,} | " + " | ".join(cells) + " |")
    lines += ["", f"Largest difference between a reported and the exact similarity: {max_error:.2e}", ""]
    return lines


def main():
    config = load_config(sys.argv[1])
    lines = [
        "# Agreement with an exhaustive exact search",
        "",
        "Each cell: queries whose hits are identical to the exhaustive reference / queries searched.",
        "Queries with hits: queries for which the reference returns at least one compound.",
        "Deviations, if any, are listed in brackets (see compare_searches.py for definitions).",
        "",
    ]
    for database in config["databases"]:
        lines += table(config, database)
    summary = "\n".join(lines)
    (EXPERIMENT_DIR / config["summary"]).write_text(summary)
    print(summary)


if __name__ == "__main__":
    main()
