"""Figure 1 — Major Land Resource Area (MLRA) 106, Nebraska portion.

Reproducible, publication-quality map of the study domain for the CSP3 LTAR
calibration site (Mead, NE), showing the Nebraska portion of MLRA 106
("Nebraska and Kansas Loess-Drift Hills", LRR M), county boundaries, the
three latitudinal subregions (Northern / Central / Southern) used to
stratify the analysis, and the dominant SSURGO soil associations within
each subregion.

Data
----
* MLRA 106 polygon — NRCS Major Land Resource Areas Geographic Database
  v5.2 (2022), queried from the authoritative NRCS FeatureServer.
* County and state boundaries — U.S. Census Bureau TIGER/Line 2023.
* Dominant soil associations — Soil Survey Staff (2020), SSURGO dominant
  component map units (association names after USDA-NRCS AH 296, 2022).
* CSP3 LTAR site — University of Nebraska ENREEC, Mead, NE (USDA-ARS
  Platte River–High Plains Aquifer LTAR).

Projection
----------
NAD83 / Conus Albers Equal Area (EPSG:5070). All distances and the scale
bar are computed in projected metres.
"""

from __future__ import annotations

import os
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import numpy as np
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, VPacker, HPacker, DrawingArea
from matplotlib.patches import Polygon as MplPolygon
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

CRS_PROJ = "EPSG:5070"       # NAD83 / Conus Albers Equal Area
CRS_GEO = "EPSG:4326"        # WGS84 lon/lat for graticule labelling

# CSP3 LTAR calibration site, ENREEC, Mead, NE (approximate centroid of the
# three CSP rotation fields; see Suyker & Verma 2012, Agric. For. Meteorol.).
CSP3_LONLAT = (-96.4766, 41.1651)

# Colour palette — colour-blind-safe, earth-toned, consistent value ramp.
COL_NORTH = "#C9A45C"   # Fillmore assoc.  (lighter, drier NCCPI band)
COL_CENTRE = "#7FA96A"  # Judson assoc.    (intermediate)
COL_SOUTH = "#3F6B46"   # Nodaway assoc.   (most productive floodplain)
COL_MLRA_EDGE = "#111111"
COL_COUNTY = "#8C8C8C"
COL_STATE = "#2E2E2E"
COL_BG = "#F2F4F5"
COL_OUTSIDE = "#E6E9EA"
COL_SUBREG = "#1C3F72"
COL_LTAR = "#B3261E"

SOIL_TABLE = [
    # (subregion, colour, dominant association, NCCPI range, brief description)
    ("Northern", COL_NORTH,
     "Fillmore assoc.",
     "0.45\u20130.65",
     "Mollic Albaqualfs, closed-depression loess uplands"),
    ("Central",  COL_CENTRE,
     "Judson assoc.",
     "0.60\u20130.75",
     "Cumulic Hapludolls, colluvial footslopes on loess"),
    ("Southern", COL_SOUTH,
     "Nodaway assoc.",
     "0.75\u20130.85",
     "Mollic Udifluvents, Holocene alluvium, Missouri R. tributaries"),
]

# ---------------------------------------------------------------------------
# Typography — Times-like serif with italic geographic feature labels
# ---------------------------------------------------------------------------

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["DejaVu Serif", "Liberation Serif", "Times New Roman",
                   "Nimbus Roman", "Times", "serif"],
    "mathtext.fontset": "dejavuserif",
    "axes.titlesize": 11,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
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
    mlra106 = gpd.read_file(DATA / "mlra_106.geojson")

    ne_counties = counties[counties["STATEFP"] == "31"].copy()
    ne_state = states[states["STUSPS"] == "NE"].copy()
    conus = states[~states["STUSPS"].isin(
        ["AK", "HI", "PR", "VI", "GU", "MP", "AS"]
    )].copy()

    for gdf in (ne_counties, ne_state, conus, mlra106):
        gdf.to_crs(CRS_PROJ, inplace=True)

    return {
        "counties": ne_counties,
        "state": ne_state,
        "conus": conus,
        "mlra106": mlra106,
    }


def build_subregions(mlra_ne: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Split the Nebraska portion of MLRA 106 into three equal-latitude
    bands (Northern / Central / Southern). The split is performed in
    geographic coordinates so the bands track true parallels."""
    mlra_geo = mlra_ne.to_crs(CRS_GEO)
    minx, miny, maxx, maxy = mlra_geo.total_bounds
    y1 = miny + (maxy - miny) / 3.0
    y2 = miny + 2.0 * (maxy - miny) / 3.0
    bands = {
        "Southern": box(minx - 1, miny - 1, maxx + 1, y1),
        "Central":  box(minx - 1, y1,       maxx + 1, y2),
        "Northern": box(minx - 1, y2,       maxx + 1, maxy + 1),
    }
    records = []
    union_mlra = unary_union(mlra_geo.geometry.values)
    for name, bnd in bands.items():
        geom = union_mlra.intersection(bnd)
        records.append({"subregion": name, "y_split_deg": (y1, y2), "geometry": geom})
    out = gpd.GeoDataFrame(records, crs=CRS_GEO).to_crs(CRS_PROJ)
    return out, (y1, y2)


# ---------------------------------------------------------------------------
# Graticule helpers
# ---------------------------------------------------------------------------

def add_graticule(ax, bounds_proj, lon_step=1.0, lat_step=0.5):
    """Draw and label parallels and meridians on a projected map."""
    minx, miny, maxx, maxy = bounds_proj
    # Generate a dense set of geographic lines, then project.
    lon_min, lat_min, lon_max, lat_max = gpd.GeoSeries(
        [box(minx, miny, maxx, maxy)], crs=CRS_PROJ
    ).to_crs(CRS_GEO).total_bounds

    lon_ticks = np.arange(np.floor(lon_min), np.ceil(lon_max) + lon_step, lon_step)
    lat_ticks = np.arange(np.floor(lat_min * 2) / 2,
                          np.ceil(lat_max * 2) / 2 + lat_step, lat_step)

    meridians = []
    for lon in lon_ticks:
        meridians.append(LineString([(lon, lat_min - 1), (lon, lat_max + 1)]))
    parallels = []
    for lat in lat_ticks:
        parallels.append(LineString([(lon_min - 1, lat), (lon_max + 1, lat)]))

    mer = gpd.GeoSeries(meridians, crs=CRS_GEO).to_crs(CRS_PROJ)
    par = gpd.GeoSeries(parallels, crs=CRS_GEO).to_crs(CRS_PROJ)
    mer.plot(ax=ax, color="#B8BEC4", linewidth=0.35, zorder=1)
    par.plot(ax=ax, color="#B8BEC4", linewidth=0.35, zorder=1)

    # Tick labels along the neatline edges.
    ax.set_xticks([])
    ax.set_yticks([])

    pad = (maxy - miny) * 0.008
    for lon, line in zip(lon_ticks, mer.geometry):
        pts = list(line.coords)
        x0, _ = pts[0]
        # Use the intersection of this meridian with the bottom of the plot.
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        # Linearly interpolate x at y = miny.
        if min(ys) <= miny <= max(ys):
            # meridians are nearly vertical in Albers; take mean x.
            x_bottom = np.interp(miny, ys, xs)
            if minx <= x_bottom <= maxx:
                ax.text(x_bottom, miny - pad, f"{abs(lon):.0f}\u00B0W",
                        ha="center", va="top", fontsize=7.5,
                        color="#333", zorder=5)

    pad_x = (maxx - minx) * 0.006
    for lat, line in zip(lat_ticks, par.geometry):
        pts = list(line.coords)
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        if min(xs) <= minx <= max(xs):
            y_left = np.interp(minx, xs, ys)
            if miny <= y_left <= maxy:
                ax.text(minx - pad_x, y_left, f"{lat:.1f}\u00B0N",
                        ha="right", va="center", fontsize=7.5,
                        color="#333", zorder=5)


# ---------------------------------------------------------------------------
# Locator inset
# ---------------------------------------------------------------------------

def add_locator_inset(fig, layers, mlra_ne, mlra_full):
    ax_in = fig.add_axes([0.065, 0.71, 0.16, 0.19])
    conus = layers["conus"]
    ne = layers["state"]
    conus.plot(ax=ax_in, facecolor="#F4EEDD", edgecolor="#9AA0A6",
               linewidth=0.35)
    # Shade Nebraska to orient the reader.
    ne.plot(ax=ax_in, facecolor="#EADCB8", edgecolor="#222", linewidth=0.6)
    # Draw the full MLRA 106 polygon (NE + KS parts) and the NE-clip.
    mlra_full.plot(ax=ax_in, facecolor=COL_LTAR, edgecolor="none",
                   alpha=0.45)
    mlra_ne.plot(ax=ax_in, facecolor=COL_LTAR, edgecolor="#5A0F0A",
                 linewidth=0.3, alpha=0.95)
    # Constrain to central CONUS for clarity.
    minx, miny, maxx, maxy = ne.total_bounds
    pad_x = (maxx - minx) * 5.5
    pad_y = (maxy - miny) * 2.0
    ax_in.set_xlim(minx - pad_x * 0.55, maxx + pad_x * 0.35)
    ax_in.set_ylim(miny - pad_y * 0.6, maxy + pad_y * 0.4)
    ax_in.set_axis_off()
    ax_in.set_aspect("equal")
    for spine in ax_in.spines.values():
        spine.set_visible(False)
    bb = ax_in.get_position()
    fig.add_artist(
        Rectangle((bb.x0, bb.y0), bb.width, bb.height,
                  transform=fig.transFigure, fill=False,
                  edgecolor="#333", linewidth=0.7, zorder=10)
    )
    ax_in.text(0.03, 0.97, "Locator",
               transform=ax_in.transAxes,
               fontsize=7.5, va="top", ha="left",
               fontstyle="italic", color="#333")
    ax_in.text(0.03, 0.05,
               "Red: MLRA 106\n(dark: NE portion)",
               transform=ax_in.transAxes,
               fontsize=6.8, va="bottom", ha="left",
               color="#333", linespacing=1.25)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> Path:
    layers = load_layers()
    ne_state = layers["state"]
    ne_counties = layers["counties"]
    mlra106 = layers["mlra106"]

    # Clip MLRA 106 to Nebraska.
    mlra_ne = gpd.overlay(mlra106, ne_state[["geometry"]], how="intersection")
    # Build subregions.
    subregions, (y1_deg, y2_deg) = build_subregions(mlra_ne)

    # Figure layout.
    fig = plt.figure(figsize=(12.0, 9.0), dpi=150, facecolor="white")
    # Main map axes.
    ax = fig.add_axes([0.055, 0.08, 0.66, 0.84])
    ax.set_facecolor(COL_BG)

    # Extent based on MLRA-NE with a generous margin.
    minx, miny, maxx, maxy = mlra_ne.total_bounds
    dx, dy = maxx - minx, maxy - miny
    pad_x, pad_y = 0.35 * dx, 0.22 * dy
    xlim = (minx - pad_x, maxx + pad_x)
    ylim = (miny - pad_y, maxy + pad_y * 1.15)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect("equal")

    # Background: all Nebraska counties within the viewport.
    # Out-of-MLRA counties are rendered in a muted neutral.
    ne_counties.plot(ax=ax, facecolor=COL_OUTSIDE,
                     edgecolor=COL_COUNTY, linewidth=0.35, zorder=2)

    # Subregions (thematic fills).
    order = ["Northern", "Central", "Southern"]
    colour_map = {"Northern": COL_NORTH, "Central": COL_CENTRE, "Southern": COL_SOUTH}
    for name in order:
        sub = subregions[subregions["subregion"] == name]
        sub.plot(ax=ax, facecolor=colour_map[name], edgecolor="none",
                 alpha=0.85, zorder=3)

    # County boundaries over fills (for geographic reference inside MLRA).
    ne_counties.boundary.plot(ax=ax, color=COL_COUNTY,
                              linewidth=0.35, zorder=4)

    # MLRA 106 outline.
    mlra_ne.boundary.plot(ax=ax, color=COL_MLRA_EDGE,
                          linewidth=1.4, zorder=6)

    # State boundary.
    ne_state.boundary.plot(ax=ax, color=COL_STATE,
                           linewidth=0.9, zorder=5)

    # Latitudinal subregion boundaries (dashed).
    for lat in (y1_deg, y2_deg):
        seg = gpd.GeoSeries(
            [LineString([(-105, lat), (-94, lat)])], crs=CRS_GEO
        ).to_crs(CRS_PROJ)
        # Clip to viewport.
        seg = seg.clip(box(*xlim, *ylim))
        seg.plot(ax=ax, color=COL_SUBREG, linewidth=1.0,
                 linestyle=(0, (6, 3)), zorder=7)

    # Graticule and tick labels.
    add_graticule(ax, (*xlim, *ylim)[:4] if False else (xlim[0], ylim[0], xlim[1], ylim[1]))

    # County name labels (major Nebraska counties overlapping the domain).
    label_counties = [
        "Dodge", "Douglas", "Sarpy", "Saunders", "Cass", "Lancaster",
        "Otoe", "Johnson", "Nemaha", "Richardson", "Pawnee", "Gage",
        "Saline", "Seward", "Jefferson",
    ]
    for _, row in ne_counties.iterrows():
        if row["NAME"] in label_counties:
            c = row.geometry.representative_point()
            ax.text(c.x, c.y, row["NAME"],
                    fontsize=7.5, style="italic", color="#333",
                    ha="center", va="center", zorder=8)

    # Subregion labels. For Northern, displace westward to avoid the
    # CSP3 LTAR star; for Central and Southern, use the interior
    # representative point.
    nudge = {
        "Northern": (-35_000, -25_000),
        "Central":  (0, 0),
        "Southern": (0, 0),
    }
    for name in order:
        sub = subregions[subregions["subregion"] == name]
        if sub.empty or sub.geometry.iloc[0].is_empty:
            continue
        c = sub.geometry.iloc[0].representative_point()
        dx_l, dy_l = nudge[name]
        ax.text(c.x + dx_l, c.y + dy_l, name,
                fontsize=11, fontweight="bold", color="#111",
                ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.28",
                          facecolor="white", edgecolor="#222",
                          alpha=0.90, linewidth=0.6),
                zorder=9)

    # CSP3 LTAR site marker.
    site = gpd.GeoSeries(
        [Point(*CSP3_LONLAT)], crs=CRS_GEO
    ).to_crs(CRS_PROJ).iloc[0]
    ax.plot(site.x, site.y, marker="*", markersize=16,
            markerfacecolor=COL_LTAR, markeredgecolor="black",
            markeredgewidth=0.6, zorder=11)
    ax.annotate(
        "CSP3 LTAR\n(Mead, NE)",
        xy=(site.x, site.y),
        xytext=(site.x + 90_000, site.y + 20_000),
        fontsize=8.5, fontweight="bold",
        color="#111",
        ha="center",
        bbox=dict(boxstyle="round,pad=0.30",
                  facecolor="white", edgecolor=COL_LTAR, linewidth=0.9),
        arrowprops=dict(arrowstyle="-", color=COL_LTAR, linewidth=0.9),
        zorder=12,
    )

    # Neatline (inner border).
    for spine in ax.spines.values():
        spine.set_color("#222")
        spine.set_linewidth(0.9)

    # Scale bar (true distance in projected metres).
    scalebar = ScaleBar(
        1.0, units="m", location="lower right",
        length_fraction=0.25, scale_loc="bottom",
        box_alpha=0.92, border_pad=0.6, pad=0.5,
        font_properties={"size": 8.5},
        color="#111", box_color="white",
        fixed_value=50, fixed_units="km",
    )
    ax.add_artist(scalebar)

    # North arrow.
    ax_na_x = xlim[1] - 0.05 * (xlim[1] - xlim[0])
    ax_na_y = ylim[0] + 0.18 * (ylim[1] - ylim[0])
    arrow_len = 0.055 * (ylim[1] - ylim[0])
    ax.annotate(
        "", xy=(ax_na_x, ax_na_y + arrow_len),
        xytext=(ax_na_x, ax_na_y),
        arrowprops=dict(arrowstyle="-|>", color="#111", linewidth=1.2,
                        mutation_scale=14),
        zorder=11,
    )
    ax.text(ax_na_x, ax_na_y + arrow_len * 1.15, "N",
            ha="center", va="bottom",
            fontsize=10, fontweight="bold", color="#111", zorder=11)

    # Legend.
    legend_handles = []
    for name, col, assoc, nccpi, _desc in SOIL_TABLE:
        legend_handles.append(
            mpatches.Patch(
                facecolor=col, edgecolor="#333", linewidth=0.4,
                label=f"{name}  \u2014  {assoc}   (NCCPI \u2248 {nccpi})",
            )
        )
    legend_handles.append(
        mlines.Line2D([], [], color=COL_SUBREG, linestyle=(0, (6, 3)),
                      linewidth=1.2, label="Latitudinal subregion boundary")
    )
    legend_handles.append(
        mlines.Line2D([], [], color=COL_MLRA_EDGE, linewidth=1.4,
                      label="MLRA 106 boundary (Nebraska portion)")
    )
    legend_handles.append(
        mpatches.Patch(facecolor=COL_OUTSIDE, edgecolor=COL_COUNTY,
                       linewidth=0.4, label="Counties outside MLRA 106")
    )
    legend_handles.append(
        mlines.Line2D([], [], color="none",
                      marker="*", markersize=10,
                      markerfacecolor=COL_LTAR, markeredgecolor="black",
                      label="CSP3 LTAR calibration site")
    )
    leg = ax.legend(
        handles=legend_handles,
        loc="lower left",
        bbox_to_anchor=(0.005, 0.005),
        title="Map legend",
        frameon=True, framealpha=0.95,
        edgecolor="#333", facecolor="white",
        borderpad=0.9, handlelength=2.0, handletextpad=0.7,
        labelspacing=0.55,
    )
    leg.get_title().set_fontweight("bold")
    leg.get_frame().set_linewidth(0.7)

    # Figure title and subtitle (single string avoids text overlap).
    fig.text(
        0.06, 0.958,
        "Figure 1.  Major Land Resource Area (MLRA) 106 \u2014 Nebraska portion",
        fontsize=13, va="top", ha="left", fontweight="bold",
    )
    fig.text(
        0.06, 0.928,
        "Study-area extent, county boundaries, latitudinal subregions, "
        "CSP3 LTAR calibration site, and dominant NCCPI soil associations",
        fontsize=10.5, va="top", ha="left", color="#333",
    )

    # Side panel: projection, datum, and data sources.
    panel_x = 0.735
    fig.text(panel_x, 0.91, "Projection & datum",
             fontsize=9.5, va="top", fontweight="bold")
    fig.text(panel_x, 0.885,
             "NAD83 / Conus Albers\nEqual Area (EPSG:5070)\n"
             "Standard parallels 29.5\u00B0N, 45.5\u00B0N\n"
             "Latitude of origin 23\u00B0N",
             fontsize=8.5, va="top", color="#222", linespacing=1.35)

    fig.text(panel_x, 0.78, "Data sources",
             fontsize=9.5, va="top", fontweight="bold")
    fig.text(panel_x, 0.755,
             "\u2022 NRCS MLRA v5.2 (2022)\n"
             "\u2022 TIGER/Line counties (2023)\n"
             "\u2022 SSURGO dominant-component\n"
             "    map units (Soil Survey Staff, 2020)\n"
             "\u2022 USDA-ARS LTAR network (CSP3)\n"
             "\u2022 National Commodity Crop\n"
             "    Productivity Index (NCCPI v3.0)",
             fontsize=8.5, va="top", color="#222", linespacing=1.35)

    fig.text(panel_x, 0.60, "Soil associations",
             fontsize=9.5, va="top", fontweight="bold")
    y = 0.575
    wrap = {
        "Fillmore assoc.":
            "Mollic Albaqualfs,\nclosed-depression\nloess uplands",
        "Judson assoc.":
            "Cumulic Hapludolls,\ncolluvial footslopes\non loess uplands",
        "Nodaway assoc.":
            "Mollic Udifluvents,\nHolocene alluvium of\nMissouri R. tributaries",
    }
    for name, col, assoc, nccpi, _desc in SOIL_TABLE:
        fig.text(panel_x, y,
                 f"{name} \u2014 {assoc}\n"
                 f"NCCPI \u2248 {nccpi}\n"
                 f"{wrap[assoc]}",
                 fontsize=8, va="top", color="#222", linespacing=1.3)
        y -= 0.095

    fig.text(panel_x, 0.19, "Cartography",
             fontsize=9.5, va="top", fontweight="bold")
    fig.text(panel_x, 0.165,
             "Equal-area statistics;\n"
             "scale bar in projected\n"
             "metres. Graticule at\n"
             "1\u00B0 / 0.5\u00B0 spacing.",
             fontsize=8.5, va="top", color="#222", linespacing=1.35)

    fig.text(panel_x, 0.085,
             "Prepared by the authors, 2026.\n"
             "Reproducible: scripts/figure1.py",
             fontsize=7.5, va="top", color="#555", fontstyle="italic",
             linespacing=1.3)

    # Locator inset.
    add_locator_inset(fig, layers, mlra_ne, mlra106)

    # Outer neatline for the whole figure.
    fig.add_artist(
        Rectangle((0.035, 0.04), 0.95, 0.935,
                  transform=fig.transFigure, fill=False,
                  edgecolor="#222", linewidth=0.9)
    )

    fig.savefig(FIG_OUT, dpi=600, bbox_inches=None,
                facecolor=fig.get_facecolor())
    fig.savefig(FIG_OUT.with_suffix(".pdf"), bbox_inches=None,
                facecolor=fig.get_facecolor())
    plt.close(fig)
    return FIG_OUT


if __name__ == "__main__":
    out = main()
    print(f"Wrote {out} ({out.stat().st_size / 1024:.1f} KB)")
