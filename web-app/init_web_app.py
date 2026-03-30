from __future__ import annotations

import argparse
import subprocess
import sys
import urllib.request
from pathlib import Path


VENDOR_ASSETS = {
    "jquery.min.js": "https://code.jquery.com/jquery-3.7.1.min.js",
    "bootstrap.min.css": "https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css",
    "bootstrap.bundle.min.js": "https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/js/bootstrap.bundle.min.js",
}


def install_flask() -> None:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "Flask>=3.0,<4"])


def download_assets(vendor_dir: Path, force: bool) -> None:
    vendor_dir.mkdir(parents=True, exist_ok=True)
    for file_name, url in VENDOR_ASSETS.items():
        destination = vendor_dir / file_name
        if destination.exists() and not force:
            print(f"Skipping existing asset: {destination}")
            continue

        print(f"Downloading {url} -> {destination}")
        with urllib.request.urlopen(url) as response, destination.open("wb") as handle:
            handle.write(response.read())


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Initialize the UFFPSim web app without Node.")
    parser.add_argument("--skip-install", action="store_true", help="Do not install Flask.")
    parser.add_argument("--skip-assets", action="store_true", help="Do not download jQuery/Bootstrap assets.")
    parser.add_argument("--force", action="store_true", help="Overwrite existing downloaded assets.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    web_app_dir = Path(__file__).resolve().parent
    vendor_dir = web_app_dir / "static" / "vendor"

    if not args.skip_install:
        install_flask()

    if not args.skip_assets:
        download_assets(vendor_dir, force=args.force)

    print("Web app initialization complete.")
    print(f"Start the app with: {sys.executable} {web_app_dir / 'app.py'} --db-file /path/to/database.h5 [--mode memory|disk] [--results ids_only|full]")


if __name__ == "__main__":
    main()