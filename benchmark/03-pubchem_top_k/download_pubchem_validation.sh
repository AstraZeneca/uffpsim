#!/usr/bin/env bash
set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
URL="https://zenodo.org/records/23057802/files/uffpsim_pubchem.tar.gz"
ARCHIVE="$DIR/uffpsim_pubchem.tar.gz"

echo "Downloading $URL ..."
curl -L -o "$ARCHIVE" "$URL"

echo "Extracting $ARCHIVE ..."
tar -xzf "$ARCHIVE" -C "$DIR"

echo "Done."
