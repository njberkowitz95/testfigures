"""Pull mapunit-level NCCPI v3 (overall) values for the dominant soil series
of MLRA 106 from the USDA-NRCS Soil Data Access (SDA) tabular service.

Outputs ``data/sda_nccpi.json`` and ``data/sda_nccpi.csv`` — the empirical
sample underlying Panel B of Figure 1. The query restricts components to
``majcompflag = 'Yes'`` and joins to the official NCCPI v3 interpretation
``'NCCPI - National Commodity Crop Productivity Index (Ver 3.0)'``.
"""

from __future__ import annotations

import csv
import json
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT_JSON = ROOT / "data" / "sda_nccpi.json"
OUT_CSV = ROOT / "data" / "sda_nccpi.csv"

SDA_URL = "https://SDMDataAccess.sc.egov.usda.gov/Tabular/post.rest"

# 13 MLRA 106 constituent-county Soil Survey Areas (Nebraska portion).
MLRA106_SSAS = [
    "NE023", "NE025", "NE067", "NE095", "NE097",
    "NE109", "NE127", "NE131", "NE133", "NE147",
    "NE151", "NE155", "NE159",
]

# Eight dominant soil series of MLRA 106 (the figure's Panel B).
SERIES = [
    "Marshall", "Monona", "Crete", "Wymore",
    "Butler", "Fillmore", "Judson", "Nodaway",
]


def query_sda() -> dict:
    areas = "','".join(MLRA106_SSAS)
    names = "','".join(SERIES)
    q = f"""
    SELECT l.areasymbol, mu.mukey, mu.muname, mu.muacres,
           c.compname, c.comppct_r, c.drainagecl, c.taxsubgrp,
           mt.interphr AS nccpi
    FROM legend l
    INNER JOIN mapunit   mu ON mu.lkey  = l.lkey
    INNER JOIN component c  ON c.mukey  = mu.mukey AND c.majcompflag = 'Yes'
    INNER JOIN cointerp  mt ON mt.cokey = c.cokey
    WHERE l.areasymbol IN ('{areas}')
      AND c.compname   IN ('{names}')
      AND mt.rulename  =
          'NCCPI - National Commodity Crop Productivity Index (Ver 3.0)'
      AND mt.mrulename =
          'NCCPI - National Commodity Crop Productivity Index (Ver 3.0)'
    """
    body = {"query": q, "format": "JSON+COLUMNNAME"}
    req = urllib.request.Request(
        SDA_URL,
        data=json.dumps(body).encode(),
        headers={"Content-Type": "application/json"},
    )
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.loads(r.read())


def main() -> None:
    data = query_sda()
    rows = data["Table"]
    cols = rows[0]
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps({"columns": cols, "rows": rows[1:]}))
    with OUT_CSV.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        w.writerows(rows[1:])
    print(f"Saved {len(rows) - 1} records \u2192 {OUT_JSON}")


if __name__ == "__main__":
    main()
