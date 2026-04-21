"""Figure 1 (revised, integrated) — Major Land Resource Area (MLRA) 106,
Nebraska portion.

Three-panel, publication-quality figure:

* **Panel A — Reference map.** The three latitudinal subregions
  (Northern / Central / Southern) shown *strictly clipped to the
  authoritative NRCS MLRA 106 polygon* (v5.2, 2022), with Nebraska county
  boundaries, a bold MLRA 106 neatline, graticule, scale bar, north
  arrow, a statewide locator inset, and the CSP3 / US-Ne3 LTAR rainfed
  calibration site at Mead, Nebraska.

* **Panel B — Empirical NCCPI v3 distributions.** Horizontal box-and-
  whisker plots of major-component NCCPI v3 *overall* index values for
  the eight dominant MLRA 106 soil series, queried live from the USDA-
  NRCS Soil Data Access (SDA) tabular service and restricted to the 13
  constituent-county Soil Survey Areas.

* **Panel C — Schematic toposequence.** Stylised hill-slope profile of
  the MLRA 106 loess-mantled till plain, with landscape-position
  assignments (summit ridgetop \u2192 floodplain) and the median NCCPI
  annotated for each series (Lewis, Pollard & Rhoades, 1967;
  USDA-NRCS, 2022).
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import numpy as np
from matplotlib.patches import Rectangle, FancyArrowPatch, Polygon as MplPolygon
from matplotlib.path import Path as MplPath
from matplotlib_scalebar.scalebar import ScaleBar
from shapely.geometry import box, LineString, Point
from shapely.ops import unary_union

# ---------------------------------------------------------------------------
# Paths and constants
# ---------------------------------------------------------------------------

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FIG_OUT = ROOT / "figures" / "Figure1.png"
FIG_OUT.parent.mkdir(parents=True, exist_ok=True)

CRS_PROJ = "EPSG:5070"   # NAD83 / Conus Albers Equal Area
CRS_GEO = "EPSG:4326"

CSP3_LONLAT = (-96.4766, 41.1651)   # US-Ne3 (rainfed), Mead, NE (ENREEC)

# Thirteen MLRA 106 constituent-county SSAs (Nebraska portion).
MLRA106_COUNTIES_NE = {
    "Saunders", "Lancaster", "Cass", "Otoe", "Johnson", "Nemaha",
    "Richardson", "Pawnee", "Gage", "Jefferson", "Saline", "Seward",
    "Butler",
}

# Latitudinal subregion palette (earth-toned, monotone value ramp).
COL_NORTH = "#E3D6B3"
COL_CENTRE = "#D9BF8D"
COL_SOUTH = "#BE9E5C"
COL_COUNTY = "#7C7C7C"
COL_COUNTY_OUT = "#C8D1D9"
COL_STATE = "#222222"
COL_MLRA_EDGE = "#111111"
COL_SUBREG = "#1F3B73"
COL_LTAR = "#A8281E"

# Drainage-class grouping (Panel B) — colour-blind-safe, matched to soil
# hydrologic behaviour rather than hue of palette in Panel A.
DRAIN_GROUPS = {
    "A": ("Well-drained loess uplands",
          ["Marshall", "Monona"], "#D9A566"),
    "B": ("Moderately well-drained uplands / footslopes",
          ["Crete", "Wymore"], "#7FA96A"),
    "C": ("Somewhat poorly drained depressions",
          ["Butler", "Fillmore"], "#B2513F"),
    "D": ("Alluvial footslopes / floodplain",
          ["Judson", "Nodaway"], "#5B90B3"),
}

SERIES_TAXO = {
    "Marshall": "Typic Hapludoll",
    "Monona":   "Typic Hapludoll",
    "Crete":    "Pachic Udertic Argiustoll",
    "Wymore":   "Aquertic Argiudoll",
    "Butler":   "Vertic Argiaquoll",
    "Fillmore": "Vertic Argialboll",
    "Judson":   "Cumulic Hapludoll",
    "Nodaway":  "Mollic Udifluvent",
}
SERIES_ORDER = ["Marshall", "Monona", "Crete", "Wymore",
                "Butler", "Fillmore", "Judson", "Nodaway"]

# NCCPI class breaks after USDA-NRCS NCCPI v3 technical note.
NCCPI_CLASSES = [
    (0.000, 0.200, "Very low", "#EEE2CC"),
    (0.200, 0.400, "Low",      "#E5D4B1"),
    (0.400, 0.600, "Mod.\nlow", "#DCC89A"),
    (0.600, 0.800, "Mod.\nhigh", "#CFB47E"),
    (0.800, 1.000, "High",     "#BFA05E"),
]

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["DejaVu Serif", "Liberation Serif", "Times New Roman",
                   "Nimbus Roman", "Times", "serif"],
    "mathtext.fontset": "dejavuserif",
    "axes.titlesize": 11,
    "axes.labelsize": 9.5,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "legend.fontsize": 8.5,
    "legend.title_fontsize": 9,
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_layers() -> dict[str, gpd.GeoDataFrame]:
    counties = gpd.read_file(DATA / "tl_2023_us_county" / "tl_2023_us_county.shp")
    states = gpd.read_file(DATA / "tl_2023_us_state" / "tl_2023_us_state.shp")
    mlra = gpd.read_file(DATA / "mlra_106.geojson")
    ne_counties = counties[counties["STATEFP"] == "31"].copy()
    ne_state = states[states["STUSPS"] == "NE"].copy()
    for gdf in (ne_counties, ne_state, mlra):
        gdf.to_crs(CRS_PROJ, inplace=True)
    return {"counties": ne_counties, "state": ne_state, "mlra": mlra}


def build_mlra_ne(mlra: gpd.GeoDataFrame, ne: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    return gpd.overlay(mlra, ne[["geometry"]], how="intersection")


def build_subregions(mlra_ne: gpd.GeoDataFrame) -> tuple[gpd.GeoDataFrame,
                                                          tuple[float, float]]:
    """Partition the NE portion of MLRA 106 into three equal-latitude bands.
    Returned bands are the *intersection of the band with the MLRA polygon*
    (not with counties), so every coloured sliver is within MLRA 106."""
    geo = mlra_ne.to_crs(CRS_GEO)
    minx, miny, maxx, maxy = geo.total_bounds
    y1 = miny + (maxy - miny) / 3.0
    y2 = miny + 2.0 * (maxy - miny) / 3.0
    bands = {
        "Southern": box(minx - 1, miny - 1, maxx + 1, y1),
        "Central":  box(minx - 1, y1,       maxx + 1, y2),
        "Northern": box(minx - 1, y2,       maxx + 1, maxy + 1),
    }
    mlra_union = unary_union(geo.geometry.values)
    out = gpd.GeoDataFrame(
        [{"subregion": name, "geometry": mlra_union.intersection(b)}
         for name, b in bands.items()],
        crs=CRS_GEO,
    ).to_crs(CRS_PROJ)
    return out, (y1, y2)


def load_sda() -> dict[str, list[float]]:
    """Return NCCPI v3 overall values per series, keyed by series name."""
    j = json.loads((DATA / "sda_nccpi.json").read_text())
    cols = j["columns"]
    out: dict[str, list[tuple[float, float]]] = defaultdict(list)
    for row in j["rows"]:
        rec = dict(zip(cols, row))
        try:
            v = float(rec["nccpi"])
            a = float(rec["muacres"] or 0.0)
        except (TypeError, ValueError):
            continue
        out[rec["compname"]].append((v, a))
    return out


# ---------------------------------------------------------------------------
# Graticule (Panel A)
# ---------------------------------------------------------------------------

def add_graticule(ax, xlim, ylim, lon_step=0.5, lat_step=0.5):
    minx, miny, maxx, maxy = xlim[0], ylim[0], xlim[1], ylim[1]
    lon_min, lat_min, lon_max, lat_max = gpd.GeoSeries(
        [box(minx, miny, maxx, maxy)], crs=CRS_PROJ
    ).to_crs(CRS_GEO).total_bounds
    lon_ticks = np.arange(np.floor(lon_min / lon_step) * lon_step,
                          np.ceil(lon_max / lon_step) * lon_step + lon_step, lon_step)
    lat_ticks = np.arange(np.floor(lat_min / lat_step) * lat_step,
                          np.ceil(lat_max / lat_step) * lat_step + lat_step, lat_step)
    meridians = gpd.GeoSeries(
        [LineString([(lon, lat_min - 1), (lon, lat_max + 1)]) for lon in lon_ticks],
        crs=CRS_GEO).to_crs(CRS_PROJ)
    parallels = gpd.GeoSeries(
        [LineString([(lon_min - 1, lat), (lon_max + 1, lat)]) for lat in lat_ticks],
        crs=CRS_GEO).to_crs(CRS_PROJ)
    meridians.plot(ax=ax, color="#BEC4CB", linewidth=0.35, zorder=1)
    parallels.plot(ax=ax, color="#BEC4CB", linewidth=0.35, zorder=1)
    ax.set_xticks([])
    ax.set_yticks([])
    # Meridian labels along the BOTTOM edge (outside map) and parallel
    # labels along the LEFT edge (outside map), to keep the interior of
    # the map clear.
    pad_y = (maxy - miny) * 0.010
    pad_x = (maxx - minx) * 0.010
    for lon, line in zip(lon_ticks, meridians.geometry):
        xs = [p[0] for p in line.coords]
        ys = [p[1] for p in line.coords]
        if min(ys) <= miny <= max(ys):
            xb = np.interp(miny, ys, xs)
            if minx + pad_x <= xb <= maxx - pad_x:
                ax.text(xb, miny - pad_y, f"{abs(lon):g}\u00B0W",
                        ha="center", va="top", fontsize=7.2,
                        color="#333", clip_on=False)
    for lat, line in zip(lat_ticks, parallels.geometry):
        xs = [p[0] for p in line.coords]
        ys = [p[1] for p in line.coords]
        if min(xs) <= minx <= max(xs):
            yl = np.interp(minx, xs, ys)
            if miny + pad_y <= yl <= maxy - pad_y:
                ax.text(minx - pad_x, yl, f"{lat:g}\u00B0N",
                        ha="right", va="center", fontsize=7.2,
                        color="#333", clip_on=False)


# ---------------------------------------------------------------------------
# Locator inset (statewide context)
# ---------------------------------------------------------------------------

def add_locator(fig, bbox, layers, mlra_ne):
    ax = fig.add_axes(bbox)
    ne = layers["state"]
    counties = layers["counties"]
    counties.plot(ax=ax, facecolor="#F6F0DF", edgecolor="#BFBFBF",
                  linewidth=0.25)
    ne.boundary.plot(ax=ax, color="#222", linewidth=0.7)
    mlra_ne.plot(ax=ax, facecolor=COL_LTAR, edgecolor="#4A0D08",
                 linewidth=0.35, alpha=0.9)
    minx, miny, maxx, maxy = ne.total_bounds
    pad_x = (maxx - minx) * 0.03
    pad_y = (maxy - miny) * 0.03
    ax.set_xlim(minx - pad_x, maxx + pad_x)
    ax.set_ylim(miny - pad_y, maxy + pad_y)
    ax.set_aspect("equal")
    ax.set_axis_off()
    for spine in ax.spines.values():
        spine.set_visible(False)
    bb = ax.get_position()
    fig.add_artist(Rectangle((bb.x0, bb.y0), bb.width, bb.height,
                             transform=fig.transFigure, fill=False,
                             edgecolor="#444", linewidth=0.7, zorder=10))
    ax.text(0.03, 0.97, "Nebraska \u2014 statewide context",
            transform=ax.transAxes, ha="left", va="top",
            fontsize=7.8, fontstyle="italic", color="#333")
    ax.text(
        0.03, 0.12,
        "MLRA 106\n(NE portion)",
        transform=ax.transAxes, ha="left", va="bottom",
        fontsize=7.5, fontweight="bold", color="#4A0D08",
        bbox=dict(boxstyle="round,pad=0.22", facecolor="white",
                  edgecolor=COL_LTAR, linewidth=0.6, alpha=0.92),
    )


# ---------------------------------------------------------------------------
# Panel A — reference map
# ---------------------------------------------------------------------------

def draw_panel_a(fig, bbox, layers, mlra_ne, subregions, split_lats):
    ax = fig.add_axes(bbox)
    ax.set_facecolor("#F3F5F6")

    ne_counties = layers["counties"]
    ne_state = layers["state"]

    # Study-area viewport around MLRA-NE with a small margin.
    minx, miny, maxx, maxy = mlra_ne.total_bounds
    dx, dy = maxx - minx, maxy - miny
    xlim = (minx - 0.30 * dx, maxx + 0.12 * dx)
    ylim = (miny - 0.12 * dy, maxy + 0.12 * dy)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect("equal")

    # County base (all Nebraska).
    ne_counties.plot(ax=ax, facecolor=COL_COUNTY_OUT,
                     edgecolor=COL_COUNTY, linewidth=0.35, zorder=2)

    # Thematic fills -- *strictly clipped to the MLRA 106 polygon*.
    colour_map = {"Northern": COL_NORTH, "Central": COL_CENTRE,
                  "Southern": COL_SOUTH}
    for name in ("Northern", "Central", "Southern"):
        band = subregions[subregions["subregion"] == name]
        if band.empty:
            continue
        band.plot(ax=ax, facecolor=colour_map[name],
                  edgecolor="none", alpha=0.95, zorder=3)

    # Re-overlay county boundaries inside MLRA for reference.
    mlra_counties = ne_counties[ne_counties["NAME"].isin(MLRA106_COUNTIES_NE)]
    mlra_counties.boundary.plot(ax=ax, color=COL_COUNTY,
                                linewidth=0.45, zorder=4)

    # MLRA 106 boundary (thick, black).
    mlra_ne.boundary.plot(ax=ax, color=COL_MLRA_EDGE,
                          linewidth=1.6, zorder=6)

    # State outline.
    ne_state.boundary.plot(ax=ax, color=COL_STATE,
                           linewidth=0.9, zorder=5)

    # Latitudinal subregion boundaries (dashed blue), clipped to viewport.
    viewport = box(xlim[0], ylim[0], xlim[1], ylim[1])
    for lat in split_lats:
        seg = gpd.GeoSeries(
            [LineString([(-105, lat), (-94, lat)])], crs=CRS_GEO
        ).to_crs(CRS_PROJ).clip(viewport)
        seg.plot(ax=ax, color=COL_SUBREG, linewidth=1.1,
                 linestyle=(0, (6, 3)), zorder=7)

    add_graticule(ax, xlim, ylim)

    # Italic county labels -- only for counties inside MLRA 106.
    for _, row in mlra_counties.iterrows():
        p = row.geometry.representative_point()
        ax.text(p.x, p.y, row["NAME"], fontsize=7.3, style="italic",
                color="#333", ha="center", va="center", zorder=8)

    # Subregion labels (displaced from CSP3 for Northern).
    nudge = {"Northern": (-40_000, -50_000),
             "Central":  (30_000, 0),
             "Southern": (10_000, 0)}
    for name in ("Northern", "Central", "Southern"):
        band = subregions[subregions["subregion"] == name]
        if band.empty or band.geometry.iloc[0].is_empty:
            continue
        c = band.geometry.iloc[0].representative_point()
        dxl, dyl = nudge[name]
        ax.text(c.x + dxl, c.y + dyl, name,
                fontsize=10.5, fontweight="bold", color="#0E1A33",
                ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.28", facecolor="white",
                          edgecolor=COL_SUBREG, linewidth=0.8, alpha=0.94),
                zorder=9)

    # CSP3 / US-Ne3 site marker (Mead, NE; rainfed CSP).
    site = gpd.GeoSeries([Point(*CSP3_LONLAT)], crs=CRS_GEO
                         ).to_crs(CRS_PROJ).iloc[0]
    ax.plot(site.x, site.y, marker="*", markersize=17,
            markerfacecolor=COL_LTAR, markeredgecolor="black",
            markeredgewidth=0.7, zorder=11)
    ax.annotate(
        "CSP3 / US-Ne3\n(Mead, NE \u2014 rainfed)",
        xy=(site.x, site.y),
        xytext=(site.x + 70_000, site.y + 15_000),
        fontsize=8.0, fontweight="bold", color="#111",
        ha="center",
        bbox=dict(boxstyle="round,pad=0.28", facecolor="white",
                  edgecolor=COL_LTAR, linewidth=0.8),
        arrowprops=dict(arrowstyle="-", color=COL_LTAR, linewidth=0.9),
        zorder=12,
    )

    # North arrow and scale bar grouped bottom-centre (beneath the map).
    x_na = xlim[0] + 0.08 * (xlim[1] - xlim[0])
    y_na = ylim[0] + 0.09 * (ylim[1] - ylim[0])
    arrow_h = 0.045 * (ylim[1] - ylim[0])
    ax.annotate("", xy=(x_na, y_na + arrow_h), xytext=(x_na, y_na),
                arrowprops=dict(arrowstyle="-|>", color="#111",
                                linewidth=1.2, mutation_scale=14),
                zorder=11)
    ax.text(x_na, y_na + arrow_h * 1.15, "N",
            ha="center", va="bottom",
            fontsize=10, fontweight="bold", color="#111", zorder=11)

    ax.add_artist(ScaleBar(
        1.0, units="m", location="lower right",
        length_fraction=0.26, scale_loc="bottom",
        box_alpha=0.94, border_pad=0.7, pad=0.5,
        font_properties={"size": 8.2},
        color="#111", box_color="white",
        fixed_value=40, fixed_units="km",
    ))

    # Neatline.
    for spine in ax.spines.values():
        spine.set_color("#222")
        spine.set_linewidth(0.9)

    # Panel letter "A" -- top-left, outside the map fills.
    ax.text(0.020, 0.975, "A", transform=ax.transAxes,
            ha="left", va="top", fontsize=13, fontweight="bold",
            color="#111",
            bbox=dict(boxstyle="square,pad=0.22", facecolor="white",
                      edgecolor="#222", linewidth=0.8))

    return ax


# ---------------------------------------------------------------------------
# Panel A legend (left of map, boxed)
# ---------------------------------------------------------------------------

def draw_panel_a_legend(fig, bbox):
    ax = fig.add_axes(bbox)
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    # Frame.
    ax.add_patch(Rectangle((0.02, 0.02), 0.96, 0.96,
                           facecolor="white", edgecolor="#333",
                           linewidth=0.8))
    ax.text(0.5, 0.93, "Map legend (Panel A)",
            ha="center", va="top", fontsize=9.2, fontweight="bold")

    entries = [
        ("swatch", COL_NORTH, "Northern subregion"),
        ("swatch", COL_CENTRE, "Central subregion"),
        ("swatch", COL_SOUTH, "Southern subregion"),
        ("dashed", COL_SUBREG, "Subregional boundary\n(N/C ; C/S)"),
        ("thick_line", COL_MLRA_EDGE, "MLRA 106 extent\n(NE portion)"),
        ("swatch", COL_COUNTY_OUT, "Counties outside\nMLRA 106"),
        ("star", COL_LTAR, "CSP3 / US-Ne3 LTAR site\n(Mead, NE \u2014 rainfed)"),
    ]
    y0 = 0.84
    dy = 0.115
    for i, (kind, col, label) in enumerate(entries):
        y = y0 - i * dy
        cx, cy = 0.12, y
        if kind == "swatch":
            ax.add_patch(Rectangle((0.08, y - 0.025), 0.10, 0.05,
                                   facecolor=col, edgecolor="#333",
                                   linewidth=0.4))
        elif kind == "dashed":
            ax.plot([0.07, 0.19], [y, y], color=col, linewidth=1.2,
                    linestyle=(0, (6, 3)))
        elif kind == "thick_line":
            ax.plot([0.07, 0.19], [y, y], color=col, linewidth=1.8)
        elif kind == "star":
            ax.plot([0.13], [y], marker="*", markersize=11,
                    markerfacecolor=col, markeredgecolor="black",
                    markeredgewidth=0.6)
        ax.text(0.24, y, label, ha="left", va="center",
                fontsize=7.8, color="#222", linespacing=1.2)


# ---------------------------------------------------------------------------
# Panel B — empirical NCCPI v3 distributions
# ---------------------------------------------------------------------------

def series_to_group(s: str) -> str:
    for g, (_desc, names, _c) in DRAIN_GROUPS.items():
        if s in names:
            return g
    return "?"


def draw_panel_b(fig, bbox, sda):
    ax = fig.add_axes(bbox)

    # Coloured NCCPI class bands (background).
    for lo, hi, _label, col in NCCPI_CLASSES:
        ax.axvspan(lo, hi, color=col, alpha=0.30, zorder=0)

    # Top-axis class-label strip.
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim() if ax.get_xlim()[1] > 0 else (0, 1))
    ax.set_xlim(0, 1)
    ax2.set_xlim(0, 1)
    for lo, hi, label, _col in NCCPI_CLASSES:
        ax2.text((lo + hi) / 2, 1.0, label,
                 transform=ax2.get_xaxis_transform(),
                 ha="center", va="bottom",
                 fontsize=7.8, fontstyle="italic", color="#444",
                 linespacing=0.9)
    ax2.set_xticks([])
    ax2.xaxis.set_visible(False)
    for spine in ax2.spines.values():
        spine.set_visible(False)

    # Box positions: top \u2192 bottom in SERIES_ORDER order.
    positions = list(range(len(SERIES_ORDER), 0, -1))
    data = []
    counts = []
    for s in SERIES_ORDER:
        vals = [v for v, _ in sda.get(s, [])]
        data.append(vals)
        acres = sum(a for _, a in sda.get(s, []))
        km2 = acres * 0.0040468564
        counts.append((len(vals), km2))

    # Colour the boxes by drainage-class group.
    colors = [DRAIN_GROUPS[series_to_group(s)][2] for s in SERIES_ORDER]

    bp = ax.boxplot(
        data, positions=positions, widths=0.52, vert=False,
        patch_artist=True, showfliers=False,
        medianprops=dict(color="#111", linewidth=1.2),
        whiskerprops=dict(color="#333", linewidth=0.9),
        capprops=dict(color="#333", linewidth=0.9),
        boxprops=dict(edgecolor="#222", linewidth=0.7),
    )
    for patch, col in zip(bp["boxes"], colors):
        patch.set_facecolor(col)
        patch.set_alpha(0.88)

    # Individual-record jitter scatter (transparency).
    rng = np.random.default_rng(42)
    for pos, vals, col in zip(positions, data, colors):
        if not vals:
            continue
        jitter = pos + rng.uniform(-0.18, 0.18, size=len(vals))
        ax.scatter(vals, jitter, s=6, color="#333", alpha=0.22,
                   linewidth=0, zorder=3)

    # Y-tick labels: series, taxonomic subgroup, and sample metadata.
    yticklabels = []
    for s in SERIES_ORDER:
        yticklabels.append(f"{s}\n({SERIES_TAXO[s]})")
    ax.set_yticks(positions)
    ax.set_yticklabels(yticklabels, fontsize=7.8, linespacing=1.2)
    ax.tick_params(axis="y", length=0)

    # Sample-size / area annotation inside the plot, to the left of each
    # box-whisker range.
    for pos, (n, km2) in zip(positions, counts):
        ax.text(0.01, pos + 0.32,
                f"$n=${n}; $A=${km2:,.0f} km\u00B2",
                ha="left", va="center",
                fontsize=6.8, color="#555", style="italic",
                zorder=4)

    ax.set_xlim(0, 1)
    ax.set_ylim(0.35, len(SERIES_ORDER) + 0.75)
    ax.set_xticks(np.arange(0, 1.01, 0.1))
    ax.set_xlabel("NCCPI v3 overall index \u2014 corn\u2013soybean productivity (0\u20131)")
    ax.grid(axis="x", linewidth=0.35, color="#CCC", zorder=0.5)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    # (Panel letter + title are drawn in figure coords by the caller.)

    # Drainage-class group legend on the right.
    handles = [
        mpatches.Patch(facecolor=DRAIN_GROUPS[g][2], edgecolor="#222",
                       linewidth=0.5,
                       label=f"Group {g} \u2014 {DRAIN_GROUPS[g][0]}")
        for g in ("A", "B", "C", "D")
    ]
    # Drainage-class labels (wrapped to fit below the axes).
    drain_labels = {
        "A": "Group A \u2014 Well-drained\n    loess uplands",
        "B": "Group B \u2014 Mod. well-drained\n    uplands/footslopes",
        "C": "Group C \u2014 Somewhat poorly\n    drained depressions",
        "D": "Group D \u2014 Alluvial footslope\n    / floodplain",
    }
    handles_wrapped = [
        mpatches.Patch(
            facecolor=DRAIN_GROUPS[g][2], edgecolor="#222", linewidth=0.5,
            label=drain_labels[g],
        )
        for g in ("A", "B", "C", "D")
    ]
    leg = ax.legend(
        handles=handles_wrapped, loc="upper center",
        bbox_to_anchor=(0.48, -0.16), ncol=4,
        title="Drainage-class group", frameon=True,
        framealpha=0.95, edgecolor="#333", borderpad=0.55,
        handlelength=1.2, labelspacing=0.35,
        columnspacing=0.7, fontsize=7.0,
        title_fontsize=8.3,
    )
    leg.get_title().set_fontweight("bold")


# ---------------------------------------------------------------------------
# Panel C — schematic toposequence
# ---------------------------------------------------------------------------

def draw_panel_c(fig, bbox, sda):
    ax = fig.add_axes(bbox)
    ax.set_xlim(0, 10)
    ax.set_ylim(-0.1, 1.15)
    ax.set_aspect("auto")
    ax.set_axis_off()
    for spine in ax.spines.values():
        spine.set_visible(False)

    # (Panel letter + title are drawn in figure coords by the caller.)

    # Medians per series (from SDA).
    med = {}
    for s in SERIES_ORDER:
        vals = [v for v, _ in sda.get(s, [])]
        med[s] = float(np.median(vals)) if vals else float("nan")

    # Terrain profile (x \u2208 [0,10]).
    xs = np.linspace(0, 10, 600)

    def terrain(x):
        y = np.piecewise(
            x,
            [x < 1.6, (x >= 1.6) & (x < 4.2),
             (x >= 4.2) & (x < 5.8), (x >= 5.8) & (x < 7.6),
             x >= 7.6],
            [
                lambda x: 0.95 - 0.02 * (x - 0.8),
                lambda x: 0.95 - 0.32 * ((x - 1.6) / 2.6),
                lambda x: 0.63 - 0.45 * ((x - 4.2) / 1.6),
                lambda x: 0.18 + 0.06 * np.sin((x - 5.8) * 2.2),
                lambda x: 0.18 - 0.10 * ((x - 7.6) / 2.4),
            ],
        )
        return y
    ys = terrain(xs)

    # Fill the profile (earth-tone).
    ax.fill_between(xs, -0.1, ys, color="#D6B97E", alpha=0.9,
                    edgecolor="#8C6A2F", linewidth=0.9, zorder=2)
    # Seasonal water table (dashed, curved).
    wt = np.where(xs < 5.8, 0.25, 0.08 + 0.02 * np.sin((xs - 5.8) * 1.8))
    ax.plot(xs, wt, color=COL_SUBREG, linestyle=(0, (5, 3)),
            linewidth=1.0, zorder=4)
    ax.text(1.1, 0.28, "Schematic seasonal water table",
            fontsize=7.5, fontstyle="italic",
            color=COL_SUBREG, va="bottom", zorder=5)

    # Wetness gradient arrow along the base.
    ax.annotate("", xy=(9.4, -0.04), xytext=(0.6, -0.04),
                arrowprops=dict(arrowstyle="->", color=COL_SUBREG,
                                linewidth=1.0))
    ax.text(5.0, -0.07, "increasing seasonal wetness",
            fontsize=7.2, fontstyle="italic",
            color=COL_SUBREG, ha="center", va="top")

    # Landscape positions \u2192 (x-centre, series, colour).
    positions = [
        (1.2, "Summit\nridgetop",
         ["Marshall", "Monona"], "A"),
        (3.2, "Interfluve\nbackslope",
         ["Crete", "Wymore"], "B"),
        (5.1, "Closed\ndepression",
         ["Butler", "Fillmore"], "C"),
        (6.7, "Footslope\n& drainageway",
         ["Judson"], "D"),
        (8.6, "Floodplain",
         ["Nodaway"], "D"),
    ]
    # NCCPI bubble vertical anchor per position (placed *above* the peak
    # for summit, slightly above the profile elsewhere).
    bubble_anchor = {
        "Summit\nridgetop":         0.55,
        "Interfluve\nbackslope":    0.44,
        "Closed\ndepression":       0.23,
        "Footslope\n& drainageway": 0.14,
        "Floodplain":               0.14,
    }
    for xc, label, series, group in positions:
        ax.text(xc, 1.02, label, ha="center", va="bottom",
                fontsize=8.5, fontweight="bold", color="#222")
        ax.text(xc, 0.96, "  \u00B7  ".join(series),
                ha="center", va="bottom",
                fontsize=7.6, fontstyle="italic", color="#333")
        y_anchor = float(terrain(np.array([xc]))[0])
        box_text = "NCCPI:\n" + "\n".join(
            f"{s} {med[s]:.2f}" for s in series
        )
        yb = bubble_anchor.get(label, y_anchor + 0.05)
        ax.text(xc, yb, box_text, ha="center", va="bottom",
                fontsize=7.2, color="#111", linespacing=1.15,
                bbox=dict(boxstyle="round,pad=0.28",
                          facecolor="white",
                          edgecolor=DRAIN_GROUPS[group][2],
                          linewidth=1.0, alpha=0.97),
                zorder=6)
        # Leader line from the box to the anchor on the profile.
        ax.plot([xc, xc], [yb, y_anchor],
                color=DRAIN_GROUPS[group][2], linewidth=0.7,
                linestyle=(0, (2, 2)), zorder=5.5)
        ax.plot([xc], [y_anchor], marker="o", markersize=5,
                markerfacecolor=DRAIN_GROUPS[group][2],
                markeredgecolor="#222", markeredgewidth=0.5,
                zorder=7)

    # y-axis label (subtle).
    ax.text(0.01, 0.55, "elevation",
            transform=ax.transAxes, rotation=90, ha="left", va="center",
            fontsize=7.5, fontstyle="italic", color="#666")


# ---------------------------------------------------------------------------
# Figure assembly
# ---------------------------------------------------------------------------

def main() -> Path:
    layers = load_layers()
    mlra_ne = build_mlra_ne(layers["mlra"], layers["state"])
    subregions, split_lats = build_subregions(mlra_ne)
    sda = load_sda()

    fig = plt.figure(figsize=(14.0, 10.0), dpi=150, facecolor="white")

    # Figure title and caption strip (top).
    fig.text(
        0.035, 0.965,
        "Figure 1.  Major Land Resource Area (MLRA) 106 \u2014 Nebraska Portion",
        fontsize=14, fontweight="bold", va="top",
    )
    fig.text(
        0.035, 0.938,
        "Geographic reference and empirical soil productivity "
        "for the rainfed maize\u2013soybean domain.",
        fontsize=10.5, va="top", color="#333", fontstyle="italic",
    )

    # ---- Layout ------------------------------------------------------------
    # Row 1: locator + Panel A + Panel B
    # Row 2: legend         + Panel C
    add_locator(fig, [0.030, 0.695, 0.120, 0.220], layers, mlra_ne)
    draw_panel_a_legend(fig, [0.020, 0.090, 0.140, 0.430])
    draw_panel_a(fig, [0.160, 0.085, 0.330, 0.830], layers, mlra_ne,
                 subregions, split_lats)
    # Panel B / C titles in figure coordinates so they don't overlap axes.
    fig.text(0.515, 0.910, "B.",
             fontsize=13, fontweight="bold", va="top", ha="left")
    fig.text(0.540, 0.910,
             "Empirical NCCPI v3 \u2014 dominant MLRA 106 series",
             fontsize=10.8, fontweight="bold", va="top", ha="left")
    draw_panel_b(fig, [0.600, 0.600, 0.320, 0.275], sda)

    fig.text(0.515, 0.475, "C.",
             fontsize=13, fontweight="bold", va="top", ha="left")
    fig.text(0.540, 0.475,
             "Schematic toposequence \u2014 MLRA 106 loess-mantled till plain",
             fontsize=10.8, fontweight="bold", va="top", ha="left")
    draw_panel_c(fig, [0.515, 0.115, 0.470, 0.335], sda)

    # Footer — data provenance block, running along the bottom.
    fig.text(
        0.035, 0.073,
        "Data sources \u2014 Panel A: NRCS MLRA v5.2 (2022) and "
        "U.S. Census TIGER/Line 2023. "
        "Panels B\u2013C: USDA-NRCS Soil Data Access (April 2026); "
        "major-component interpretations restricted to the 13\n"
        "MLRA 106 constituent-county SSAs (Saunders, Cass, Lancaster, Otoe, "
        "Johnson, Nemaha, Richardson, Pawnee, Gage, Jefferson, Saline, "
        "Seward, Butler). "
        "Projection: NAD83 / Conus Albers Equal Area (EPSG 5070).\n"
        "Landscape-position assignments after Lewis, Pollard & Rhoades "
        "(1967) and USDA-NRCS (2022); box-plot medians denote the empirical "
        "median NCCPI v3 overall value across all major components.",
        fontsize=7.3, va="top", color="#555", linespacing=1.45,
    )

    # Outer neatline.
    fig.add_artist(Rectangle((0.015, 0.015), 0.975, 0.965,
                             transform=fig.transFigure, fill=False,
                             edgecolor="#222", linewidth=0.9))

    fig.savefig(FIG_OUT, dpi=600, bbox_inches=None,
                facecolor=fig.get_facecolor())
    fig.savefig(FIG_OUT.with_suffix(".pdf"), bbox_inches=None,
                facecolor=fig.get_facecolor())
    plt.close(fig)
    return FIG_OUT


if __name__ == "__main__":
    out = main()
    print(f"Wrote {out} ({out.stat().st_size / 1024:.1f} KB)")
