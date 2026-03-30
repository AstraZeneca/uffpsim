from __future__ import annotations

import argparse
import base64
import io
import json
import threading
from pathlib import Path
from typing import Any

from flask import Flask, jsonify, render_template, request
from rdkit import Chem
from rdkit.Chem import Draw

import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
	sys.path.insert(0, str(REPO_ROOT))

from uffpsim import UFFPSimSearchEngine, get_database_info


def create_app(
	db_file: str | Path | None = None,
	search_mode: str = "memory",
	results_mode: str = "ids_only",
) -> Flask:
	app = Flask(__name__, template_folder="templates", static_folder="static")
	app.config["UFFPSIM_REPO_ROOT"] = REPO_ROOT
	app.config["DB_FILE"] = None
	app.config["DB_INFO"] = {}
	app.config["SEARCH_MODE"] = search_mode
	app.config["RESULTS_MODE"] = results_mode
	app.config["FP_TYPE"] = None
	app.config["FP_BITS_SIZE"] = None
	app.config["NUM_MOLS"] = None
	app.config["CLUSTER_THRESHOLD"] = None
	app.search_engine = None
	app.search_lock = threading.Lock()
	app.search_queue_depth = 0
	app.search_queue_lock = threading.Lock()

	if db_file is not None:
		resolved_db_file = Path(db_file).expanduser().resolve()
		if not resolved_db_file.is_file():
			raise FileNotFoundError(f"Database file not found: {resolved_db_file}")
		app.config["DB_FILE"] = str(resolved_db_file)
		app.config["DB_INFO"] = _safe_database_info(resolved_db_file)
		app.search_engine = UFFPSimSearchEngine(str(resolved_db_file), mode=search_mode)
		if results_mode == "full":
			app.search_engine.build_mol_id_to_index_map()
		_fp_meta = json.loads(app.search_engine.fp_store.fp_params_json)
		app.config["FP_TYPE"] = _fp_meta.get("fp_type")
		app.config["FP_BITS_SIZE"] = app.search_engine.fp_store.fp_bits_size
		app.config["NUM_MOLS"] = app.search_engine.num_mols
		app.config["CLUSTER_THRESHOLD"] = app.search_engine.fp_store.inner_clustering_threshold

	@app.get("/")
	def index() -> str:
		return render_template("index.html")

	@app.get("/api/config")
	def config() -> Any:
		return jsonify(
			{
				"repo_root": str(REPO_ROOT),
				"db_file": app.config["DB_FILE"],
				"db_info": app.config["DB_INFO"],
				"search_mode": app.config["SEARCH_MODE"],
				"results_mode": app.config["RESULTS_MODE"],
				"fp_type": app.config["FP_TYPE"],
				"fp_bits_size": app.config["FP_BITS_SIZE"],
				"num_mols": app.config["NUM_MOLS"],
				"cluster_threshold": app.config["CLUSTER_THRESHOLD"],
				"queue_depth": app.search_queue_depth,
			}
		)

	@app.post("/api/search")
	def search() -> Any:
		if app.search_engine is None or app.config["DB_FILE"] is None:
			return jsonify({"error": "Search engine is not initialized. Start the server with --db-file."}), 400

		try:
			payload = request.get_json(silent=True) or {}
			smiles_list = _parse_smiles_list(payload)
			if not smiles_list:
				return jsonify({"error": "Provide at least one SMILES string."}), 400

			threshold = float(payload.get("threshold", 0.6))
			limit_by = int(payload.get("limit_by", 10))
		except Exception as exc:
			return jsonify({"error": str(exc)}), 400

		# Track how many requests are waiting + running
		with app.search_queue_lock:
			app.search_queue_depth += 1
		try:
			with app.search_lock:  # only one search runs at a time; others wait
				try:
					engine = app.search_engine
					search_results = engine.batch_search(smiles_list, threshold, limit_by)
				except Exception as exc:
					return jsonify({"error": str(exc)}), 400
		finally:
			with app.search_queue_lock:
				app.search_queue_depth -= 1

		results_mode = app.config["RESULTS_MODE"]
		response = []
		for query_smiles, hits in zip(smiles_list, search_results):
			query_record = {
				"query_smiles": query_smiles,
				"query_image": _smiles_to_data_url(query_smiles),
				"hits": [],
				"error": None,
			}
			if hits is None:
				query_record["error"] = "Invalid SMILES or fingerprint generation failed."
				response.append(query_record)
				continue

			for compound_id, score in hits:
				hit: dict[str, Any] = {"compound_id": compound_id, "score": float(score)}
				if results_mode == "full":
					hit_smiles = engine.get_smiles_for_id(compound_id)
					hit["smiles"] = hit_smiles
					hit["image"] = _smiles_to_data_url(hit_smiles)
				query_record["hits"].append(hit)
			response.append(query_record)

		return jsonify(
			{
				"db_file": app.config["DB_FILE"],
				"threshold": threshold,
				"limit_by": limit_by,
				"results": response,
			}
		)

	return app


def _safe_database_info(db_file: Path) -> dict[str, Any]:
	try:
		return get_database_info(str(db_file))
	except Exception:
		return {}


def _parse_smiles_list(payload: dict[str, Any]) -> list[str]:
	smiles_list = payload.get("smiles_list")
	if isinstance(smiles_list, list):
		return [str(item).strip() for item in smiles_list if str(item).strip()]

	smiles_text = str(payload.get("smiles_text") or "")
	return [line.strip() for line in smiles_text.splitlines() if line.strip()]


def _smiles_to_data_url(smiles: str) -> str | None:
	mol = Chem.MolFromSmiles(smiles)
	if mol is None:
		return None

	image = Draw.MolToImage(mol, size=(280, 180))
	buffer = io.BytesIO()
	image.save(buffer, format="PNG")
	encoded = base64.b64encode(buffer.getvalue()).decode("ascii")
	return f"data:image/png;base64,{encoded}"


app = create_app()


def parse_args() -> argparse.Namespace:
	parser = argparse.ArgumentParser(description="Run the UFFPSim web app with a preloaded database.")
	parser.add_argument("--db-file", required=True, help="Path to the uffpsim .h5 database file to load at startup.")
	parser.add_argument("--mode", default="memory", choices=["memory", "disk"], help="Search engine load mode (default: memory).")
	parser.add_argument(
		"--results",
		default="ids_only",
		choices=["ids_only", "full"],
		help="Result detail level: 'ids_only' returns compound IDs and scores only; "
			 "'full' also returns SMILES strings and RDKit images (builds an ID→SMILES map at startup). "
			 "Default: ids_only.",
	)
	parser.add_argument("--host", default="127.0.0.1", help="Host interface to bind.")
	parser.add_argument("--port", type=int, default=5000, help="TCP port to bind.")
	parser.add_argument("--debug", action="store_true", help="Enable Flask debug mode.")
	return parser.parse_args()


if __name__ == "__main__":
	args = parse_args()
	app = create_app(args.db_file, search_mode=args.mode, results_mode=args.results)
	app.run(host=args.host, port=args.port, debug=args.debug)
