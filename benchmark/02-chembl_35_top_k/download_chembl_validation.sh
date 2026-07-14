#!/usr/bin/env bash
set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
URL="https://zenodo.org/records/21171543/files/uffpsim_chembl_35.tar.gz"
ARCHIVE="$DIR/uffpsim_chembl_35.tar.gz"

echo "Downloading $URL ..."
curl -L -o "$ARCHIVE" "$URL"

echo "Extracting $ARCHIVE ..."
tar -xzf "$ARCHIVE" -C "$DIR"

echo "Done."
