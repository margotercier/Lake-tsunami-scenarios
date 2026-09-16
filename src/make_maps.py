"""
Hazard maps and figures for the Lake Hawea landslide-tsunami scenarios.

Colour use follows a strict sequential scheme: one hue, light to dark, for each
magnitude field (blue for wave height on the lake, orange for inundation depth on
land - two sequential contexts, so the second takes the next hue). No rainbow ramps.
"""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, LightSource
import matplotlib.patheffects as pe

DATA, OUT = "/home/user/Lake-tsunami-scenarios/data", "/home/user/Lake-tsunami-scenarios/outputs"
FIG = "/home/user/Lake-tsunami-scenarios/outputs/figures"
os.makedirs(FIG, exist_ok=True)
RES, LAKE_LEVEL, XMIN, YMAX = 30.0, 345.0, 1292000.0, 5096000.0
SUB = (80, 1520, 100, 990)

COARSEN = 5          # must match seiche_modes.py
BLUE  = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
ORANGE = ["#fde3d6", "#f9bfa4", "#f39a74", "#eb6834", "#c94f1f", "#a03d16", "#782d0f"]
ORDINAL3 = ["#86b6ef", "#2a78d6", "#104281"]     # ordinal blue steps 250/450/650
CMAP_W = LinearSegmentedColormap.from_list("wave", BLUE)
CMAP_L = LinearSegmentedColormap.from_list("inund", ORANGE)
INK, INK2 = "#0b0b0b", "#52514e"

PLACES = {"Lake Hāwea township": (1302581, 5053435), "Hāwea dam": (1302244, 5052629)}
LABEL_OFF = {"Lake Hāwea township": (8, 6), "Hāwea dam": (-10, -14)}
LABEL_HA  = {"Lake Hāwea township": "left", "Hāwea dam": "right"}
SCENARIOS = ["S1_moderate", "S2_large", "S3_extreme"]
LABELS = {"S1_moderate": "S1  moderate  5×10⁵ m³",
          "S2_large":    "S2  large  5×10⁶ m³",
          "S3_extreme":  "S3  extreme  2×10⁷ m³"}


def sub(a):
    r0, r1, c0, c1 = SUB
    return a[r0:r1, c0:c1]


def to_rc(x, y):
    return (YMAX - y) / RES - SUB[0], (x - XMIN) / RES - SUB[2]


def load():
    dem = sub(np.load(f"{DATA}/dem.npy"))
    lake = sub(np.load(f"{DATA}/lake_mask.npy"))
    depth = sub(np.load(f"{DATA}/depth.npy"))
    return dem, lake, depth


def hillshade(dem):
    ls = LightSource(azdeg=315, altdeg=45)
    return ls.hillshade(dem, vert_exag=1.5, dx=RES, dy=RES)


def base_axes(ax, hs, lake, extent):
    ax.imshow(hs, cmap="gray", vmin=0, vmax=1.35, extent=extent, origin="upper",
              interpolation="bilinear")
    ax.imshow(np.where(lake, 1.0, np.nan), cmap=LinearSegmentedColormap.from_list(
        "w", ["#e8eef5", "#e8eef5"]), extent=extent, origin="upper", alpha=0.55)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_color("#c9c8c3")


def annotate_places(ax, fs=8, which=None):
    for name, (x, y) in PLACES.items():
        if which and name not in which:
            continue
        r, c = to_rc(x, y)
        ax.plot(c, r, "o", ms=5, mfc="#ffffff", mec=INK, mew=1.4, zorder=6)
        ax.annotate(name, (c, r), textcoords="offset points",
                    xytext=LABEL_OFF.get(name, (7, 4)),
                    ha=LABEL_HA.get(name, "left"),
                    fontsize=fs, color=INK, zorder=6,
                    path_effects=[pe.withStroke(linewidth=2.5, foreground="white")])


# ------------------------------------------------------------------ figure 1
def fig_overview(dem, lake, depth, hs):
    ny, nx = dem.shape
    ext = (0, nx, ny, 0)
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 9.2), constrained_layout=True)
    vmax = 0.0
    fields = {}
    for s in SCENARIOS:
        me = np.load(f"{OUT}/{s}_maxeta.npy")
        wave = np.where(lake & np.isfinite(me), me - LAKE_LEVEL, np.nan)
        zb = np.where(lake, LAKE_LEVEL - depth, dem)
        inund = np.where((~lake) & np.isfinite(me) & (me > zb + 0.1), me - zb, np.nan)
        fields[s] = (wave, inund)
        vmax = max(vmax, np.nanpercentile(wave, 99.5))
    vmax = float(np.ceil(vmax))

    for ax, s in zip(axes, SCENARIOS):
        wave, inund = fields[s]
        base_axes(ax, hs, lake, ext)
        im = ax.imshow(wave, cmap=CMAP_W, norm=Normalize(0, vmax), extent=ext,
                       origin="upper", interpolation="nearest")
        il = ax.imshow(inund, cmap=CMAP_L, norm=Normalize(0, 8), extent=ext,
                       origin="upper", interpolation="nearest")
        meta = json.load(open(f"{OUT}/{s}_meta.json"))
        r, c = to_rc(*meta["impact_xy"])
        ax.plot(c, r, marker="*", ms=17, mfc="#e34948", mec="white", mew=1.2, zorder=7)
        annotate_places(ax)
        ax.set_title(LABELS[s], fontsize=11, color=INK, pad=8)
        ax.text(0.03, 0.015, f"max wave {np.nanmax(wave):.0f} m\n"
                             f"inundated {np.isfinite(inund).sum()*RES*RES/1e6:.2f} km²",
                transform=ax.transAxes, fontsize=8.5, color=INK2, va="bottom",
                bbox=dict(fc="white", ec="none", alpha=.8, pad=3))
    cb2 = fig.colorbar(il, ax=axes, location="bottom", fraction=.035, pad=.02, aspect=55)
    cb2.set_label("maximum inundation depth on land (m)", fontsize=9)
    cb1 = fig.colorbar(im, ax=axes, location="bottom", fraction=.035, pad=.01, aspect=55)
    cb1.set_label("maximum wave height on the lake (m above normal level)", fontsize=9)
    for cb in (cb1, cb2):
        cb.ax.tick_params(labelsize=8, colors=INK2)
    fig.get_layout_engine().set(rect=(0, 0, 1, 0.935))
    fig.suptitle("Lake Hāwea landslide-tsunami scenarios — modelled maximum wave field",
                 fontsize=14, color=INK, y=0.995)
    fig.text(0.5, 0.952, "Alpine Fault–triggered rock avalanche entering the lake at the "
             "★ (mid-lake east shore, the highest-ranked source zone)",
             ha="center", fontsize=9.5, color=INK2, va="top")
    fig.savefig(f"{FIG}/01_overview.png", dpi=150, facecolor="white")
    plt.close(fig)
    print("wrote 01_overview.png")


# ------------------------------------------------------------------ figure 2
def fig_township(dem, lake, depth, hs):
    r, c = to_rc(*PLACES["Lake Hāwea township"])
    pad_y, pad_x = 150, 170
    r0, r1 = int(r - pad_y), int(r + pad_y // 2)
    c0, c1 = int(c - pad_x), int(c + pad_x)
    ext = (c0, c1, r1, r0)
    fig, axes = plt.subplots(1, 3, figsize=(14, 6.4), constrained_layout=True)
    zb_full = np.where(lake, LAKE_LEVEL - depth, dem)
    for ax, s in zip(axes, SCENARIOS):
        me = np.load(f"{OUT}/{s}_maxeta.npy")
        inund = np.where((~lake) & np.isfinite(me) & (me > zb_full + 0.1), me - zb_full, np.nan)
        wave = np.where(lake & np.isfinite(me), me - LAKE_LEVEL, np.nan)
        ax.imshow(hs[r0:r1, c0:c1], cmap="gray", vmin=0, vmax=1.35, extent=ext,
                  origin="upper", interpolation="bilinear")
        ax.imshow(np.where(lake, 1.0, np.nan)[r0:r1, c0:c1],
                  cmap=LinearSegmentedColormap.from_list("w", ["#e8eef5"] * 2),
                  extent=ext, origin="upper", alpha=.55)
        ax.imshow(wave[r0:r1, c0:c1], cmap=CMAP_W, norm=Normalize(0, 12), extent=ext,
                  origin="upper")
        il = ax.imshow(inund[r0:r1, c0:c1], cmap=CMAP_L, norm=Normalize(0, 8),
                       extent=ext, origin="upper")
        a = np.isfinite(inund[r0:r1, c0:c1]).sum() * RES * RES / 1e6
        ax.set_title(f"{LABELS[s]}\ninundated here: {a:.2f} km²", fontsize=10, color=INK)
        annotate_places(ax, fs=9)
        ax.set_xticks([]); ax.set_yticks([])
        # 1 km scale bar
        ax.plot([c0 + 12, c0 + 12 + 1000 / RES], [r1 - 14] * 2, "-", lw=3, color=INK)
        ax.text(c0 + 12, r1 - 20, "1 km", fontsize=8, color=INK)
    cb = fig.colorbar(il, ax=axes, location="bottom", fraction=.045, pad=.02, aspect=50)
    cb.set_label("maximum inundation depth on land (m)", fontsize=9)
    cb.ax.tick_params(labelsize=8, colors=INK2)
    fig.suptitle("Lake Hāwea township — modelled inundation", fontsize=14, color=INK)
    fig.savefig(f"{FIG}/02_township.png", dpi=160, facecolor="white")
    plt.close(fig)
    print("wrote 02_township.png")


# ------------------------------------------------------------------ figure 3
def fig_arrival(lake, hs):
    ny, nx = lake.shape
    ext = (0, nx, ny, 0)
    arr = np.load(f"{OUT}/S2_large_arrival.npy")
    arr = np.where(lake, arr, np.nan) / 60.0
    fig, ax = plt.subplots(figsize=(6.6, 9.4), constrained_layout=True)
    base_axes(ax, hs, lake, ext)
    im = ax.imshow(arr, cmap=CMAP_W, norm=Normalize(0, 12), extent=ext, origin="upper")
    cs = ax.contour(arr, levels=[1, 2, 3, 4, 5, 6, 8, 10], colors="#0b0b0b",
                    linewidths=.7, extent=ext, origin="upper")
    ax.clabel(cs, fmt="%d min", fontsize=7.5)
    meta = json.load(open(f"{OUT}/S2_large_meta.json"))
    r, c = to_rc(*meta["impact_xy"])
    ax.plot(c, r, marker="*", ms=18, mfc="#e34948", mec="white", mew=1.2, zorder=7)
    annotate_places(ax, fs=9)
    cb = fig.colorbar(im, ax=ax, fraction=.04, pad=.02)
    cb.set_label("wave arrival time (minutes after slide impact)", fontsize=9)
    cb.ax.tick_params(labelsize=8, colors=INK2)
    ax.set_title("Wave arrival time — scenario S2 (5×10⁶ m³)", fontsize=12.5, color=INK)
    fig.savefig(f"{FIG}/03_arrival.png", dpi=150, facecolor="white")
    plt.close(fig)
    print("wrote 03_arrival.png")


# ------------------------------------------------------------------ figure 4
def fig_gauges():
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True, constrained_layout=True)
    keys = ["Hawea township", "Hawea dam"]
    titles = ["Lake Hāwea township foreshore", "Hāwea dam / lake outlet"]
    for ax, key, tt in zip(axes, keys, titles):
        for s, col in zip(SCENARIOS, ORDINAL3):
            m = json.load(open(f"{OUT}/{s}_meta.json"))
            t = np.array(m["times"]) / 60.0
            y = np.array(m["gauges"][key])
            ax.plot(t, y, lw=2, color=col, label=LABELS[s], solid_capstyle="round")
            ax.annotate(LABELS[s].split()[1], (t[-1], y[-1]), fontsize=8, color=col,
                        xytext=(4, 0), textcoords="offset points", va="center")
        ax.axhline(0, color="#c9c8c3", lw=1)
        ax.set_title(tt, fontsize=11, color=INK, loc="left")
        ax.set_ylabel("water level anomaly (m)", fontsize=9, color=INK2)
        ax.grid(True, color="#ecebe7", lw=.8)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[-1].set_xlabel("time after slide impact (minutes)", fontsize=9, color=INK2)
    axes[0].legend(frameon=False, fontsize=9, loc="upper right")
    fig.suptitle("Modelled water-level time series", fontsize=13, color=INK)
    fig.savefig(f"{FIG}/04_gauges.png", dpi=160, facecolor="white")
    plt.close(fig)
    print("wrote 04_gauges.png")


# ------------------------------------------------------------------ figure 5
def fig_seiche(lake, hs):
    info = json.load(open(f"{DATA}/seiche_modes.json"))
    ny, nx = np.load(f"{DATA}/dem.npy").shape          # full-grid shape, pre-SUB
    fig, axes = plt.subplots(1, 4, figsize=(13, 6.6), constrained_layout=True)
    div = LinearSegmentedColormap.from_list("div", ["#184f95", "#6da7ec", "#f0efec",
                                                    "#f39a74", "#a03d16"])
    for k, ax in enumerate(axes, start=1):
        # modes are solved on a 5x-coarsened grid: upsample back before slicing
        coarse = np.load(f"{DATA}/seiche_mode{k}.npy")
        full = np.kron(coarse, np.ones((COARSEN, COARSEN), dtype="float32"))
        pad = np.full((ny, nx), np.nan, dtype="float32")
        h_, w_ = min(full.shape[0], ny), min(full.shape[1], nx)
        pad[:h_, :w_] = full[:h_, :w_]
        sh = pad[SUB[0]:SUB[1], SUB[2]:SUB[3]]
        ax.imshow(hs, cmap="gray", vmin=0, vmax=1.4, origin="upper", interpolation="bilinear")
        ax.imshow(sh, cmap=div, vmin=-1, vmax=1, origin="upper", interpolation="nearest")
        T = info["modes"][k - 1]["period_min"]
        ax.set_title(f"mode {k}\nT = {T:.1f} min", fontsize=10.5, color=INK)
        ax.set_xticks([]); ax.set_yticks([])
    fig.get_layout_engine().set(rect=(0, 0.10, 1, 0.93))
    fig.suptitle("Lake Hāwea natural seiche modes (free oscillations of the basin)",
                 fontsize=13, color=INK, y=0.995)
    fig.text(.5, .035, "Blue and orange are opposite phases of the standing wave; the "
             "shoreline between them is a node.\nStrong, long-period shaking in an Alpine "
             "Fault rupture can excite these modes, which then ring for hours.",
             ha="center", fontsize=9, color=INK2)
    fig.savefig(f"{FIG}/05_seiche.png", dpi=150, facecolor="white")
    plt.close(fig)
    print("wrote 05_seiche.png")


# ------------------------------------------------------------------ figure 6
def fig_volume_height():
    """Modelled scenarios against the historical landslide-tsunami record."""
    hist = [("Deep Cove, Fiordland 1987", 0.1, 2.5), ("Gold Arm, Fiordland 2003", 0.2, 4.5),
            ("Mt Colonel Foster, Canada 1946", 0.7, 29), ("Chehalis Lake, Canada 2007", 2.5, 25),
            ("Tafjord, Norway 1934", 3.0, 62), ("Lake Askja, Iceland 2014", 20.0, 71),
            ("Lituya Bay, Alaska 1958", 30.0, 524)]
    fig, ax = plt.subplots(figsize=(9.5, 6.6), constrained_layout=True)
    v = np.array([h[1] for h in hist]); r = np.array([h[2] for h in hist])
    ax.scatter(v, r, s=64, facecolor="white", edgecolor=INK2, lw=1.5, zorder=4,
               label="observed landslide-generated waves")
    for n, vv, rr in hist:
        ax.annotate(n, (vv, rr), fontsize=7.6, color=INK2, xytext=(7, -3),
                    textcoords="offset points")
    k = np.polyfit(np.log10(v), np.log10(r), 1)
    xs = np.logspace(-2, 2, 50)
    ax.plot(xs, 10 ** np.polyval(k, np.log10(xs)), color="#c9c8c3", lw=1.6, ls="--",
            zorder=2, label=f"power-law fit  R ∝ V^{k[0]:.2f}")
    for s, col in zip(SCENARIOS, ORDINAL3):
        m = json.load(open(f"{OUT}/{s}_meta.json"))
        ax.scatter([m["scenario"]["volume"] / 1e6], [m["H_M"]], s=130, color=col,
                   zorder=6, marker="D", edgecolor="white", lw=1.4)
        ax.annotate(f'  {s.split("_")[0]}: near-field {m["H_M"]:.0f} m', 
                    (m["scenario"]["volume"] / 1e6, m["H_M"]), fontsize=9, color=col,
                    xytext=(9, 4), textcoords="offset points", weight="bold")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("landslide volume (million m³)", fontsize=10, color=INK2)
    ax.set_ylabel("maximum wave / run-up height (m)", fontsize=10, color=INK2)
    ax.grid(True, which="both", color="#ecebe7", lw=.8); ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    ax.set_title("How big is a lake tsunami? Modelled Hāwea scenarios against the "
                 "observed record", fontsize=12.5, color=INK, loc="left")
    fig.savefig(f"{FIG}/06_volume_vs_height.png", dpi=160, facecolor="white")
    plt.close(fig)
    print("wrote 06_volume_vs_height.png")


# ------------------------------------------------------------------ figure 7
def fig_tectonic(dem, lake, depth, hs):
    """Hunter Valley Fault rupture: whole-lake field and township detail."""
    if not os.path.exists(f"{OUT}/S4_tectonic_maxeta.npy"):
        print("skip 07: S4 not run"); return
    me = np.load(f"{OUT}/S4_tectonic_maxeta.npy")
    meta = json.load(open(f"{OUT}/S4_tectonic_meta.json"))
    zb = np.where(lake, LAKE_LEVEL - depth, dem)
    wave = np.where(lake & np.isfinite(me), me - LAKE_LEVEL, np.nan)
    inund = np.where((~lake) & np.isfinite(me) & (me > zb + 0.1), me - zb, np.nan)
    ny, nx = dem.shape
    fig = plt.figure(figsize=(12.5, 8.6), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.25])
    ax = fig.add_subplot(gs[0, 0])
    base_axes(ax, hs, lake, (0, nx, ny, 0))
    im = ax.imshow(wave, cmap=CMAP_W, norm=Normalize(0, max(1.0, np.nanpercentile(wave, 99.5))),
                   extent=(0, nx, ny, 0), origin="upper")
    # mapped fault trace
    from pyproj import Transformer
    fl = json.load(open(f"{DATA}/faults_near_hawea.json"))["Hunter Valley Fault"]["pts"]
    tr = Transformer.from_crs(4326, 2193, always_xy=True)
    fxs, fys = tr.transform(np.array(fl)[:, 0], np.array(fl)[:, 1])
    rr, cc = to_rc(np.array(fxs), np.array(fys))
    ax.plot(cc, rr, ".", ms=2.0, color="#e34948", zorder=7)
    ax.plot([], [], "-", color="#e34948", lw=2, label="Hunter Valley Fault (NZAFD)")
    ax.legend(frameon=False, fontsize=8.5, loc="upper left")
    annotate_places(ax, fs=8.5)
    cb = fig.colorbar(im, ax=ax, fraction=.04, pad=.02)
    cb.set_label("maximum wave height (m)", fontsize=9)
    cb.ax.tick_params(labelsize=8, colors=INK2)
    ax.set_title(f"S4 · Hunter Valley Fault rupture\n{meta['uz_max']:+.1f} m / "
                 f"{meta['uz_min']:+.1f} m lakebed displacement", fontsize=11, color=INK)

    ax2 = fig.add_subplot(gs[0, 1])
    m2 = json.load(open(f"{OUT}/S4_tectonic_meta.json"))
    t = np.array(m2["times"]) / 60.0
    for key, col, lab in (("Hawea township", ORDINAL3[1], "Lake Hāwea township"),
                          ("Hawea dam", ORDINAL3[2], "Hāwea dam")):
        y = np.array(m2["gauges"][key])
        ax2.plot(t, y, lw=2, color=col, label=lab)
        ax2.annotate(lab, (t[-1], y[-1]), fontsize=8.5, color=col, xytext=(4, 0),
                     textcoords="offset points", va="center")
    ax2.axhline(0, color="#c9c8c3", lw=1)
    ax2.set_xlabel("time after rupture (minutes)", fontsize=9, color=INK2)
    ax2.set_ylabel("water level anomaly (m)", fontsize=9, color=INK2)
    ax2.grid(True, color="#ecebe7", lw=.8); ax2.set_axisbelow(True)
    for sp in ("top", "right"):
        ax2.spines[sp].set_visible(False)
    ax2.legend(frameon=False, fontsize=9)
    ax2.set_title("Water level at the southern end", fontsize=11, color=INK, loc="left")
    fig.suptitle("Scenario S4 — co-seismic rupture of the Hunter Valley Fault beneath "
                 "Lake Hāwea", fontsize=13, color=INK)
    fig.savefig(f"{FIG}/07_tectonic.png", dpi=150, facecolor="white")
    plt.close(fig)
    print("wrote 07_tectonic.png")


if __name__ == "__main__":
    dem, lake, depth = load()
    hs = hillshade(dem)
    fig_overview(dem, lake, depth, hs)
    fig_township(dem, lake, depth, hs)
    fig_arrival(lake, hs)
    fig_gauges()
    fig_seiche(lake, hs)
    fig_volume_height()
    fig_tectonic(dem, lake, depth, hs)
