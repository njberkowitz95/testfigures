# testfigures

Reproducible, publication-quality figures for the Major Land Resource Area
(MLRA) 106 — Nebraska portion — CSP3 LTAR study.

## Figure 1

![Figure 1](figures/Figure1.png)

### Academic caption

> **Figure 1.** Spatial context of the study domain, Major Land Resource
> Area (MLRA) 106, *Nebraska and Kansas Loess–Drift Hills* (Land Resource
> Region M), restricted to the Nebraska portion. The domain is delineated
> using the NRCS *Major Land Resource Areas Geographic Database*
> v5.2 (USDA–NRCS, 2022) and clipped to the Nebraska state polygon from
> TIGER/Line 2023 (U.S. Census Bureau, 2023). County geometries are
> likewise drawn from TIGER/Line 2023. For the purpose of trend
> stratification, the domain is partitioned into three equal-latitude
> subregions (Northern, Central, Southern), whose boundaries are shown as
> dashed lines. The dominant soil associations within each subregion are
> derived from SSURGO dominant-component map units (Soil Survey Staff,
> 2020): *Fillmore association* (Mollic Albaqualfs, closed-depression
> loess uplands; NCCPI ≈ 0.45–0.65) in the Northern subregion,
> *Judson association* (Cumulic Hapludolls, colluvial footslopes on loess;
> NCCPI ≈ 0.60–0.75) in the Central subregion, and *Nodaway association*
> (Mollic Udifluvents on Holocene alluvium of Missouri-River tributaries;
> NCCPI ≈ 0.75–0.85) in the Southern subregion, where NCCPI values
> summarise the National Commodity Crop Productivity Index v3.0 (USDA–
> NRCS, 2020). The CSP3 LTAR calibration site (USDA–ARS Platte River–High
> Plains Aquifer LTAR; Mead, Nebraska, 41.17 °N, 96.48 °W) is indicated
> by the red star. All geospatial operations (area calculation, latitude
> banding, clipping, and the printed scale bar) are performed in the
> NAD83 / Conus Albers Equal Area projection (EPSG 5070). The inset in
> the upper-left locates the study area within the contiguous United
> States; the red polygon delineates MLRA 106 in its entirety, with the
> darker shade denoting the Nebraska portion used in this study.

### Reproduce

```bash
pip install -r requirements.txt
# Source data (places into ./data):
bash scripts/fetch_data.sh
# Render figure:
python3 scripts/figure1.py
```

The script writes both `figures/Figure1.png` (600 dpi) and
`figures/Figure1.pdf` (vector) for submission-grade output.

### Layers and sources

| Layer | Source | Citation |
| --- | --- | --- |
| MLRA 106 polygon | NRCS *Major Land Resource Areas* v5.2 FeatureServer | USDA–NRCS (2022) |
| County boundaries | U.S. Census Bureau TIGER/Line 2023 (`tl_2023_us_county`) | U.S. Census Bureau (2023) |
| State boundaries | U.S. Census Bureau TIGER/Line 2023 (`tl_2023_us_state`) | U.S. Census Bureau (2023) |
| Dominant soil associations | SSURGO dominant-component map units | Soil Survey Staff (2020) |
| CSP3 LTAR site | ENREEC, University of Nebraska | Suyker & Verma (2012) |

### Cartographic specification

* Projection: NAD83 / Conus Albers Equal Area (EPSG 5070)
* Graticule: 1° meridians, 0.5° parallels
* Scale bar: 50 km, computed in projected metres
* Typography: serif, 8–13 pt; italic labels for geographic features
* Palette: earth-toned, colour-blind-safe, with monotone value ramp
  aligned to the NCCPI gradient (lighter = lower, darker = higher)
