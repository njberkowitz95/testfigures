#!/usr/bin/env bash
# Fetch authoritative geospatial inputs for Figure 1.
#
# Outputs:
#   data/tl_2023_us_county/*.shp   (U.S. Census TIGER/Line 2023)
#   data/tl_2023_us_state/*.shp    (U.S. Census TIGER/Line 2023)
#   data/mlra_106.geojson          (NRCS MLRA v5.2 2022; MLRARSYM = 106)

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DATA="$ROOT/data"
mkdir -p "$DATA"

# --- TIGER/Line counties ----------------------------------------------------
if [ ! -f "$DATA/tl_2023_us_county/tl_2023_us_county.shp" ]; then
  curl -sL -o "$DATA/tl_2023_us_county.zip" \
    "https://www2.census.gov/geo/tiger/TIGER2023/COUNTY/tl_2023_us_county.zip"
  unzip -oq "$DATA/tl_2023_us_county.zip" -d "$DATA/tl_2023_us_county"
  rm -f "$DATA/tl_2023_us_county.zip"
fi

# --- TIGER/Line states ------------------------------------------------------
if [ ! -f "$DATA/tl_2023_us_state/tl_2023_us_state.shp" ]; then
  curl -sL -o "$DATA/tl_2023_us_state.zip" \
    "https://www2.census.gov/geo/tiger/TIGER2023/STATE/tl_2023_us_state.zip"
  unzip -oq "$DATA/tl_2023_us_state.zip" -d "$DATA/tl_2023_us_state"
  rm -f "$DATA/tl_2023_us_state.zip"
fi

# --- MLRA 106 polygon (NRCS authoritative FeatureServer) --------------------
if [ ! -f "$DATA/mlra_106.geojson" ]; then
  curl -sG \
    --data-urlencode "where=MLRARSYM='106'" \
    --data-urlencode "outFields=MLRA_ID,MLRARSYM,MLRA_NAME,LRRSYM" \
    --data-urlencode "outSR=4326" \
    --data-urlencode "f=geojson" \
    "https://services.arcgis.com/SXbDpmb7xQkk44JV/arcgis/rest/services/Major_Land_Resource_Areas/FeatureServer/0/query" \
    -o "$DATA/mlra_106.geojson"
fi

echo "Done. Inputs in $DATA"
