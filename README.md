# testfigures

Reproducible, publication-quality figures for the Major Land Resource Area
(MLRA) 106 — Nebraska portion — CSP3 / US-Ne3 LTAR rainfed maize–soybean
study.

## Figure 1 (integrated, three-panel)

![Figure 1](figures/Figure1.png)

Figure 1 is an integrated, three-panel figure:

* **Panel A — Reference map.** Latitudinal subregions of the study
  domain (Northern / Central / Southern) are **strictly clipped to the
  authoritative NRCS MLRA 106 polygon** (v5.2, 2022), not to county
  outlines, so the thematic fills correspond unambiguously to MLRA
  extent. Nebraska counties inside MLRA 106 are labelled in italic; the
  MLRA boundary is drawn as a bold black neatline. The CSP3 / US-Ne3
  LTAR rainfed calibration site (Mead, NE; 41.17 °N, 96.48 °W) is
  marked, and a statewide locator inset anchors the map within
  Nebraska. Projection: NAD83 / Conus Albers Equal Area (EPSG 5070).
* **Panel B — Empirical NCCPI v3 distributions.** Horizontal box-and-
  whisker plot of *major-component* NCCPI v3 overall index values for
  the eight dominant MLRA 106 soil series (Marshall, Monona, Crete,
  Wymore, Butler, Fillmore, Judson, Nodaway), pulled live from the
  USDA-NRCS Soil Data Access (SDA) tabular service and restricted to
  the 13 constituent-county Soil Survey Areas (SSAs) of MLRA 106.
  Boxes are coloured by drainage-class group (A = well-drained loess
  uplands, B = moderately well-drained uplands/footslopes, C =
  somewhat poorly drained depressions, D = alluvial footslope /
  floodplain). Individual major-component records are overlaid as
  jitter points for transparency.
* **Panel C — Schematic toposequence.** Stylised hill-slope profile of
  the MLRA 106 loess-mantled till plain. Landscape positions (summit
  ridgetop → floodplain, left to right) are annotated with the series
  typical of each position and the **empirical median NCCPI v3 value
  for each series** taken from Panel B. Catena relationships follow
  Lewis, Pollard & Rhoades (1967) and USDA-NRCS (2022).

### Academic caption

> **Figure 1.** Major Land Resource Area (MLRA) 106 — Nebraska portion:
> geographic reference and empirical soil productivity for the rainfed
> maize–soybean domain. (A) NRCS MLRA 106 polygon (v5.2, 2022), clipped
> to Nebraska, partitioned into three equal-latitude subregions
> (Northern / Central / Southern) and annotated with Nebraska counties
> (italic) and the CSP3 / US-Ne3 LTAR rainfed calibration site at Mead,
> NE. Projection: NAD83 / Conus Albers Equal Area (EPSG 5070). (B)
> Distribution of NCCPI v3 overall index values (0–1) for the eight
> dominant soil series of MLRA 106, restricted to major components of
> the 13 constituent-county SSAs (Saunders, Cass, Lancaster, Otoe,
> Johnson, Nemaha, Richardson, Pawnee, Gage, Jefferson, Saline, Seward,
> Butler), queried from the USDA-NRCS Soil Data Access service; each
> box is coloured by drainage-class group and annotated with the
> sample size *n* and total map-unit area *A*. (C) Schematic
> toposequence of the MLRA 106 loess-mantled till plain, showing the
> five dominant landscape positions (summit ridgetop → interfluve
> backslope → closed depression → footslope / drainageway → floodplain)
> and the empirical median NCCPI v3 value for the series typical of
> each position; catena relationships follow Lewis, Pollard & Rhoades
> (1967) and USDA-NRCS (2022).

### Reproduce

```bash
pip install -r requirements.txt
bash scripts/fetch_data.sh          # MLRA 106 polygon + TIGER/Line 2023
python3 scripts/fetch_sda_nccpi.py  # SDA query (live, ~30 KB)
python3 scripts/figure1.py          # renders PNG (600 dpi) + PDF
```

### Layers and sources

| Layer | Source | Citation |
| --- | --- | --- |
| MLRA 106 polygon | NRCS *Major Land Resource Areas* v5.2 FeatureServer | USDA-NRCS (2022) |
| County boundaries | U.S. Census TIGER/Line 2023 (`tl_2023_us_county`) | U.S. Census (2023) |
| State boundaries | U.S. Census TIGER/Line 2023 (`tl_2023_us_state`) | U.S. Census (2023) |
| NCCPI v3 overall | USDA-NRCS Soil Data Access (SDA), `cointerp` table | Soil Survey Staff (accessed 2026) |
| CSP3 / US-Ne3 site | UNL-ENREEC rainfed CSP rotation, Mead, NE | Suyker & Verma (2012) |

### Cartographic specification

* Projection: NAD83 / Conus Albers Equal Area (EPSG 5070)
* Graticule: 0.5° spacing on both axes
* Scale bar: 40 km, computed in projected metres
* Typography: serif, 7–14 pt; italic labels for geographic features
* Palette: monotone value ramp for latitudinal subregions (low → high
  NCCPI, light → dark), and a four-colour qualitative palette for
  drainage-class groups (Panel B / C)
