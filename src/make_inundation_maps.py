"""
"What gets covered" maps: zoomed inundation extent at the places that matter.

Each panel shows nested flood extents - the area covered in the smallest scenario is
drawn darkest, with the extra ground each larger scenario adds shown progressively
lighter. Darker therefore means "wet even in the modest scenario", i.e. more certain.
"""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import LinearSegmentedColormap, LightSource
from matplotlib.patches import Patch
from scipy import ndimage

DATA = "/home/user/Lake-tsunami-scenarios/data"
OUT  = "/home/user/Lake-tsunami-scenarios/outputs"
FIG  = f"{OUT}/figures"
RES, LAKE_LEVEL, XMIN, YMAX = 30.0, 345.0, 1292000.0, 5096000.0
SUB = (80, 1520, 100, 990)
SCEN = ["S1_moderate", "S2_large", "S3_extreme"]
# ordinal blue: darkest = flooded even in S1
BAND = {"S1_moderate": "#0d366b", "S2_large": "#2a78d6", "S3_extreme": "#9ec5f4"}
LAKEC, INK, INK2 = "#dbe7f2", "#0b0b0b", "#52514e"

AREAS = [
    dict(key="township", title="Lake Hāwea township, outlet and dam",
         xy=(1302490, 5053400), half=(1800, 1500),
         marks={"township": (1302581, 5053435)},
         marks2={"outlet / dam": (1302485, 5053490)},
         note="The township sits on a terrace above the lake; what floods is the\n"
              "foreshore, the boat ramp and the outlet channel below the dam."),
    dict(key="head", title="Head of the lake — Hunter River delta",
         xy=(1316111, 5089088), half=(2700, 2700), marks={},
         note="The flattest ground on the lake. Shallow water spreads a long way\n"
              "inland here — the largest area covered anywhere around Hāwea."),
    dict(key="slide", title="The slide zone — mid-lake east shore",
         xy=(1308080, 5069930), half=(2600, 2200),
         marks={"slide impact": (1308080, 5069930)},
         note="Run-up of tens of metres onto steep ground. Little flat land to\n"
              "cover, but nothing on this shore would survive it."),

]


def load():
    sub = lambda a: a[SUB[0]:SUB[1], SUB[2]:SUB[3]]
    dem = sub(np.load(f"{DATA}/dem.npy"))
    lake = sub(np.load(f"{DATA}/lake_mask.npy"))
    depth = sub(np.load(f"{DATA}/depth.npy"))
    zbl = sub(np.load(f"{DATA}/zb_land.npy"))
    zb = np.where(lake, LAKE_LEVEL - depth, zbl)
    return dem, lake, zb


def rc(x, y):
    return (YMAX - y) / RES - SUB[0], (x - XMIN) / RES - SUB[2]


def main():
    dem, lake, zb = load()
    ls = LightSource(azdeg=315, altdeg=45)
    hs = ls.hillshade(dem, vert_exag=1.6, dx=RES, dy=RES)

    floods = {}
    for s in SCEN:
        me = np.load(f"{OUT}/{s}_maxeta.npy")
        floods[s] = (~lake) & np.isfinite(me) & (me > zb + 0.1)

    stats = {}
    fig, axes = plt.subplots(2, 2, figsize=(13.6, 12.4), constrained_layout=True)
    panels = list(axes.ravel())
    overview_ax = panels.pop(2)          # bottom-left becomes the whole-lake view
    for ax, area in zip(panels, AREAS):
        r, c = rc(*area["xy"])
        hy, hx = area["half"][1] / RES, area["half"][0] / RES
        r0, r1 = int(max(0, r - hy)), int(min(dem.shape[0], r + hy))
        c0, c1 = int(max(0, c - hx)), int(min(dem.shape[1], c + hx))
        ext = (c0, c1, r1, r0)
        sl = (slice(r0, r1), slice(c0, c1))

        ax.imshow(hs[sl], cmap="gray", vmin=0.05, vmax=1.35, extent=ext,
                  origin="upper", interpolation="bilinear")
        ax.imshow(np.where(lake[sl], 1.0, np.nan), extent=ext, origin="upper",
                  cmap=LinearSegmentedColormap.from_list("l", [LAKEC, LAKEC]))
        # nested bands, largest first so the smallest ends up on top
        for s in reversed(SCEN):
            m = floods[s][sl]
            ax.imshow(np.where(m, 1.0, np.nan), extent=ext, origin="upper",
                      cmap=LinearSegmentedColormap.from_list("b", [BAND[s], BAND[s]]),
                      interpolation="nearest")
        # crisp shoreline
        ax.contour(lake[sl].astype(float), levels=[0.5], colors="#5b6570",
                   linewidths=.7, extent=ext, origin="upper")

        for nm, (x, y) in list(area["marks"].items()) + list(area.get("marks2", {}).items()):
            below = nm in area.get("marks2", {})
            mr, mc = rc(x, y)
            ax.plot(mc, mr, "o" if "slide" not in nm else "*",
                    ms=7 if "slide" not in nm else 17,
                    mfc="white" if "slide" not in nm else "#e34948",
                    mec=INK if "slide" not in nm else "white", mew=1.4, zorder=8)
            ax.annotate(nm, (mc, mr), xytext=(-9, -14) if below else (8, 5),
                        ha="right" if below else "left",
                        textcoords="offset points", fontsize=9, color=INK, zorder=8,
                        path_effects=[pe.withStroke(linewidth=2.6, foreground="white")])

        a3 = floods["S3_extreme"][sl].sum() * RES * RES / 1e6
        a1 = floods["S1_moderate"][sl].sum() * RES * RES / 1e6
        a2 = floods["S2_large"][sl].sum() * RES * RES / 1e6
        stats[area["key"]] = dict(title=area["title"], S1=round(a1, 2),
                                  S2=round(a2, 2), S3=round(a3, 2))
        s1txt = "not reached in S1" if a1 < 0.005 else f"{a1:.2f} km² in S1"
        ax.set_title(f"{area['title']}\n{s1txt}  ·  {a3:.2f} km² covered in S3",
                     fontsize=10.5, color=INK)
        ax.text(0.03, 0.03, area["note"], transform=ax.transAxes, fontsize=8.5,
                color=INK2, va="bottom", linespacing=1.5,
                bbox=dict(fc="white", ec="#dcdfdb", alpha=.92, boxstyle="round,pad=0.45"))
        # 1 km scale bar
        sbx = c1 - 10 - 1000 / RES
        ax.plot([sbx, sbx + 1000 / RES], [r1 - 14] * 2, "-", lw=3.2, color=INK,
                solid_capstyle="butt",
                path_effects=[pe.withStroke(linewidth=5, foreground="white")])
        ax.text(sbx + 500 / RES, r1 - 22, "1 km", fontsize=8.5, color=INK, ha="center",
                path_effects=[pe.withStroke(linewidth=2.5, foreground="white")])
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_color("#c9c8c3")

    # ---- whole-lake overview, with the zoom boxes marked -------------------
    ax = overview_ax
    ny, nx = dem.shape
    ext = (0, nx, ny, 0)
    ax.imshow(hs, cmap="gray", vmin=0.05, vmax=1.35, extent=ext, origin="upper",
              interpolation="bilinear")
    ax.imshow(np.where(lake, 1.0, np.nan), extent=ext, origin="upper",
              cmap=LinearSegmentedColormap.from_list("l", [LAKEC, LAKEC]))
    for s_ in reversed(SCEN):
        # dilate so a one-cell strip is still visible at this zoom
        m = ndimage.binary_dilation(floods[s_], np.ones((3, 3)))
        ax.imshow(np.where(m, 1.0, np.nan), extent=ext, origin="upper",
                  cmap=LinearSegmentedColormap.from_list("b", [BAND[s_], BAND[s_]]),
                  interpolation="nearest")
    for i, area in enumerate(AREAS, start=1):
        r, c = rc(*area["xy"])
        hy, hx = area["half"][1] / RES, area["half"][0] / RES
        ax.add_patch(plt.Rectangle((c - hx, r - hy), 2 * hx, 2 * hy, fill=False,
                                   ec="#a8391b", lw=1.6, zorder=9))
        ax.annotate(str(i), (c - hx + 4, r - hy + 4), fontsize=11, color="white",
                    weight="bold", va="top", zorder=10,
                    bbox=dict(fc="#a8391b", ec="none", boxstyle="circle,pad=0.22"))
    tot = floods["S3_extreme"].sum() * RES * RES / 1e6
    ax.set_title(f"The whole lake — every shore that gets wet\n"
                 f"{tot:.1f} km² covered in total (S3); boxes 1-3 are the panels above",
                 fontsize=10.5, color=INK)
    ax.text(0.03, 0.975, "Most of Hāwea's shoreline is too steep for water to get\n"
            "far inland. The exceptions are the deltas — and those are\n"
            "where the flat, usable land is.", transform=ax.transAxes, fontsize=8.5,
            color=INK2, va="top", linespacing=1.5,
            bbox=dict(fc="white", ec="#dcdfdb", alpha=.92, boxstyle="round,pad=0.45"))
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_color("#c9c8c3")

    handles = [Patch(fc=BAND[s], ec="none",
                     label={"S1_moderate": "covered even in S1 (5×10⁵ m³)",
                            "S2_large": "also covered in S2 (5×10⁶ m³)",
                            "S3_extreme": "also covered in S3 (2×10⁷ m³)"}[s])
               for s in SCEN] + [Patch(fc=LAKEC, ec="#5b6570", label="lake at normal level")]
    fig.legend(handles=handles, loc="outside lower center", ncol=4, frameon=False,
               fontsize=9.5)
    fig.suptitle("Lake Hāwea — what the tsunami actually covers", fontsize=15, color=INK)
    stats["_total"] = {s_: round(floods[s_].sum() * RES * RES / 1e6, 2) for s_ in SCEN}
    json.dump(stats, open(f"{OUT}/areas_covered.json", "w"), indent=2)
    fig.savefig(f"{FIG}/08_areas_covered.png", dpi=155, facecolor="white")
    plt.close(fig)
    print("wrote 08_areas_covered.png and areas_covered.json")


if __name__ == "__main__":
    main()
