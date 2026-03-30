from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

from .database import (
    build_mol_id_index_table,
    create_database,
    create_database_parallel,
    redo_inner_clustering,
)
from .search_engine import UFFPSimSearchEngine


def _json_or_string(value: str) -> dict[str, Any] | str:
    try:
        return json.loads(value)
    except json.JSONDecodeError:
        return value


def _json_dict_or_none(value: str | None) -> dict[str, Any] | None:
    if value is None:
        return None
    parsed = json.loads(value)
    if parsed is None:
        return None
    if not isinstance(parsed, dict):
        raise argparse.ArgumentTypeError("Expected a JSON object.")
    return parsed


def _add_create_database_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser("create-database", help="Create a new uffpsim HDF5 database.")
    parser.add_argument("-i", "--input-file", required=True, help="Input molecule file path.")
    parser.add_argument("-d", "--db-file", required=True, help="Output HDF5 database path.")
    parser.add_argument("-f", "--fp-type", required=True, help="Fingerprint type.")
    parser.add_argument("-p", "--fp-params", type=_json_dict_or_none, default=None, help="Fingerprint parameters as JSON object.")
    parser.add_argument("-w", "--workers", type=int, default=1, help="Number of workers. Use >1 for parallel supplier mode.")
    parser.add_argument("-g", "--gen-ids", action=argparse.BooleanOptionalAction, default=True, help="Generate IDs if missing in input.")
    parser.add_argument("-m", "--mol-id-prop", default=None, help="Property name to use as molecule ID.")
    parser.add_argument("-l", "--mol-id-max-chars", type=int, default=15, help="Maximum molecule ID length.")
    parser.add_argument("-I", "--info", default="", help="Additional database info as JSON string or plain text.")
    parser.add_argument("-t", "--inner-clustering-threshold", type=float, default=0.2, help="Inner clustering threshold.")
    parser.add_argument("-c", "--cluster-mode", choices=["memory", "disk"], default="memory", help="Clustering mode.")
    parser.add_argument("-P", "--cluster-parallel", action="store_true", help="Enable OpenMP-based clustering parallelism.")
    parser.set_defaults(command=_run_create_database)


def _add_redo_clustering_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser("redo-clustering", help="Redo inner clustering for an existing database.")
    parser.add_argument("-d", "--db-file", required=True, help="HDF5 database path.")
    parser.add_argument("-t", "--threshold", required=True, type=float, help="New inner clustering threshold.")
    parser.add_argument("-c", "--cluster-mode", choices=["memory", "disk"], default="memory", help="Clustering mode.")
    parser.add_argument("-P", "--cluster-parallel", action="store_true", help="Enable OpenMP-based clustering parallelism.")
    parser.set_defaults(command=_run_redo_clustering)


def _add_build_index_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser("build-mol-id-index-table", help="Build serialized MolIdIndexTable for an existing database.")
    parser.add_argument("-d", "--db-file", required=True, help="HDF5 database path.")
    parser.set_defaults(command=_run_build_mol_id_index_table)


def _add_launch_web_app_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser("launch-web-app", help="Launch the uffpsim web application.")
    parser.add_argument("-d", "--db-file", required=True, help="HDF5 database path to preload.")
    parser.add_argument("-m", "--mode", choices=["memory", "disk"], default="memory", help="Search engine load mode.")
    parser.add_argument("-r", "--results", choices=["ids_only", "full"], default="ids_only", help="Result detail mode.")
    parser.add_argument("-H", "--host", default="127.0.0.1", help="Host interface to bind.")
    parser.add_argument("-p", "--port", type=int, default=5000, help="TCP port to bind.")
    parser.add_argument("-D", "--debug", action="store_true", help="Enable Flask debug mode.")
    parser.set_defaults(command=_run_launch_web_app)


def _add_search_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser("search", help="Search one or more SMILES and write hits to CSV.")
    parser.add_argument("-d", "--db-file", required=True, help="HDF5 database path.")
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("-q", "--smiles", help="Single SMILES query from command line.")
    input_group.add_argument("-i", "--smiles-file", help="Path to input file containing one SMILES per line.")
    parser.add_argument("-o", "--output-csv", required=True, help="Output CSV file path.")
    parser.add_argument("-m", "--mode", choices=["memory", "disk"], default="memory", help="Search engine load mode.")
    parser.add_argument("-t", "--threshold", type=float, default=0.6, help="Similarity threshold.")
    parser.add_argument("-k", "--limit", type=int, default=10, help="Maximum hits per query.")
    parser.add_argument("-s", "--include-hit-smiles", action="store_true", help="Include hit SMILES in output CSV.")
    parser.set_defaults(command=_run_search_to_csv)


def _load_smiles_queries(args: argparse.Namespace) -> list[str]:
    if args.smiles is not None:
        smiles = args.smiles.strip()
        if not smiles:
            raise ValueError("--smiles cannot be empty")
        return [smiles]

    if args.smiles_file is None:
        raise ValueError("Provide either --smiles or --smiles-file")

    smiles_path = Path(args.smiles_file)
    if not smiles_path.is_file():
        raise FileNotFoundError(f"SMILES input file not found: {smiles_path}")

    smiles_list: list[str] = []
    with smiles_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            candidate = line.strip()
            if not candidate or candidate.startswith("#"):
                continue
            smiles_list.append(candidate)

    if not smiles_list:
        raise ValueError("No SMILES found in input file")
    return smiles_list


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="uffpsim", description="UFFPSim command-line interface.")
    subparsers = parser.add_subparsers(dest="subcommand", required=True)

    _add_create_database_parser(subparsers)
    _add_redo_clustering_parser(subparsers)
    _add_build_index_parser(subparsers)
    _add_search_parser(subparsers)
    _add_launch_web_app_parser(subparsers)
    return parser


def _run_create_database(args: argparse.Namespace) -> int:
    info = _json_or_string(args.info)
    common_kwargs = {
        "input_file": args.input_file,
        "db_file": args.db_file,
        "fp_type": args.fp_type,
        "fp_params": args.fp_params,
        "gen_ids": args.gen_ids,
        "mol_id_prop": args.mol_id_prop,
        "mol_id_max_chars": args.mol_id_max_chars,
        "info": info,
        "inner_clustering_threshold": args.inner_clustering_threshold,
        "cluster_mode": args.cluster_mode,
        "cluster_parallel": args.cluster_parallel,
    }
    if args.workers > 1:
        create_database_parallel(workers=args.workers, **common_kwargs)
    else:
        create_database(**common_kwargs)
    return 0


def _run_redo_clustering(args: argparse.Namespace) -> int:
    redo_inner_clustering(
        db_file=args.db_file,
        threshold=args.threshold,
        cluster_mode=args.cluster_mode,
        cluster_parallel=args.cluster_parallel,
    )
    return 0


def _run_build_mol_id_index_table(args: argparse.Namespace) -> int:
    build_mol_id_index_table(args.db_file)
    return 0


def _run_search_to_csv(args: argparse.Namespace) -> int:
    smiles_queries = _load_smiles_queries(args)

    engine = UFFPSimSearchEngine(args.db_file, mode=args.mode)
    if args.include_hit_smiles:
        engine.build_mol_id_to_index_map()

    results = engine.batch_search(smiles_queries, args.threshold, args.limit)

    output_path = Path(args.output_csv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["query_index", "query_smiles", "hit_rank", "compound_id", "score", "hit_smiles", "error"])

        for query_index, (query_smiles, hits) in enumerate(zip(smiles_queries, results), start=1):
            if hits is None:
                writer.writerow([query_index, query_smiles, "", "", "", "", "Invalid SMILES or fingerprint generation failed"])
                continue

            if len(hits) == 0:
                writer.writerow([query_index, query_smiles, "", "", "", "", "No hits"])
                continue

            for hit_rank, (compound_id, score) in enumerate(hits, start=1):
                hit_smiles = ""
                if args.include_hit_smiles:
                    hit_smiles = engine.get_smiles_for_id(compound_id)
                writer.writerow([query_index, query_smiles, hit_rank, compound_id, float(score), hit_smiles, ""])

    print(f"Wrote search results to: {output_path}")
    return 0


def _run_launch_web_app(args: argparse.Namespace) -> int:
    from .web_app.app import create_app as create_web_app

    app = create_web_app(args.db_file, search_mode=args.mode, results_mode=args.results)
    app.run(host=args.host, port=args.port, debug=args.debug)
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.command(args)


if __name__ == "__main__":
    raise SystemExit(main())