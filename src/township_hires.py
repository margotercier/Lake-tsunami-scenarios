"""
High-resolution inundation map of Lake Hawea township.

Combines three open datasets, all fetched directly over HTTP:
  * LINZ 0.1 m urban aerial imagery (Queenstown Lakes 2022-2023) as the basemap,
  * LINZ 1 m LiDAR bare-earth DEM (Wanaka 2022-2023) for ground elevation,
  * the modelled peak water level from the 30 m shallow-water model.

The flood extent is a connectivity-constrained fill: a 1 m cell is covered if the
ground is below the modelled peak water level AND is connected to the lake across
ground that is also below it. This is downscaling of a modelled level, not a 1 m
hydrodynamic simulation - it gives the extent of ground the modelled peak reaches.
"""
import os, json
import numpy as np
import rasterio
from rasterio.windows import from_bounds
from scipy import ndimage

os.environ.setdefault("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
os.environ.setdefault("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", ".tiff,.tif")

DATA = "/home/user/Lake-tsunami-scenarios/data"
OUT  = "/home/user/Lake-tsunami-scenarios/outputs"
S3   = "/vsicurl/https://nz-elevation.s3.ap-southeast-2.amazonaws.com/otago/wanaka_2022-2023/dem_1m/2193"
IMG  = "/vsicurl/https://nz-imagery.s3.ap-southeast-2.amazonaws.com/otago/queenstown-lakes_2022-2023_0.1m/rgb/2193"
LIDAR_TILES = ["CA13_10000_0401", "CA13_10000_0501"]
RES30, LAKE_LEVEL, XMIN, YMAX = 30.0, 345.0, 1292000.0, 5096000.0
SUB = (80, 1520, 100, 990)
SCEN = ["S1_moderate", "S2_large", "S3_extreme"]
# Lake frontage only. The Hawea River corridor below the outlet is deliberately out
# of scope: the LiDAR shows the channel bed ~20 m below lake level with no dam
# structure resolved in it, so anything this model said about downstream flooding
# would be an artefact of where the barrier is assumed to sit, not a result.
WIN = (1302150, 5053380, 1304400, 5053860)
MPP = 0.5                      # working resolution, metres per pixel


def mosaic(paths, win, mpp, bands=1):
    x0, y0, x1, y1 = win
    W, H = int((x1 - x0) / mpp), int((y1 - y0) / mpp)
    shape = (H, W) if bands == 1 else (bands, H, W)
    out = np.full(shape, np.nan, dtype="float32")
    for p in paths:
        try:
            ds = rasterio.open(p)
        except Exception as e:
            print(f"  skip {p.rsplit('/',1)[-1]}: {e}"); continue
        with ds:
            b = ds.bounds
            xx0, yy0 = max(x0, b.left), max(y0, b.bottom)
            xx1, yy1 = min(x1, b.right), min(y1, b.top)
            if xx1 <= xx0 or yy1 <= yy0:
                continue
            oh, ow = int(round((yy1 - yy0) / mpp)), int(round((xx1 - xx0) / mpp))
            if oh < 1 or ow < 1:
                continue
            idx = list(range(1, min(bands, ds.count) + 1))
            a = ds.read(idx, window=from_bounds(xx0, yy0, xx1, yy1, ds.transform),
                        out_shape=(len(idx), oh, ow), boundless=False).astype("float32")
            if ds.nodata is not None:
                a[a == ds.nodata] = np.nan
            a[a < -1000] = np.nan
            r0, c0 = int(round((y1 - yy1) / mpp)), int(round((xx0 - x0) / mpp))
            if bands == 1:
                sl = out[r0:r0 + oh, c0:c0 + ow]
                m = np.isfinite(a[0]); sl[m] = a[0][m]
            else:
                sl = out[:, r0:r0 + oh, c0:c0 + ow]
                m = np.isfinite(a)
                n = min(sl.shape[0], a.shape[0])
                sl[:n][m[:n]] = a[:n][m[:n]]
    return out


def upto(a, shape, k):
    """Nearest-neighbour upsample by k, then crop/pad to exactly `shape`."""
    b = np.kron(a, np.ones((k, k), dtype=a.dtype))
    out = np.full(shape, b.flat[0] if b.size else 0, dtype=b.dtype)
    h, w = min(shape[0], b.shape[0]), min(shape[1], b.shape[1])
    out[:h, :w] = b[:h, :w]
    if h < shape[0]:
        out[h:, :w] = b[b.shape[0] - 1, :w]
    if w < shape[1]:
        out[:, w:] = out[:, w - 1][:, None]
    return out


def main():
    x0, y0, x1, y1 = WIN
    print(f"window {x1-x0:.0f} x {y1-y0:.0f} m at {MPP} m/px")

    print("fetching 1 m LiDAR ...")
    dem = mosaic([f"{S3}/{t}.tiff" for t in LIDAR_TILES], WIN, MPP)
    print(f"  LiDAR valid {100*np.isfinite(dem).mean():.1f}%")

    print("fetching 0.1 m aerial imagery ...")
    tiles = json.load(open("/tmp/town_img_tiles.json"))
    rgb = mosaic([f"{IMG}/{t}.tiff" for t in tiles], WIN, MPP, bands=3)
    print(f"  imagery valid {100*np.isfinite(rgb[0]).mean():.1f}%")

    c0, c1 = int((x0 - XMIN) / RES30), int((x1 - XMIN) / RES30)
    r0, r1 = int((YMAX - y1) / RES30), int((YMAX - y0) / RES30)
    k = int(RES30 / MPP)

    # --- modelled peak water level, nearest lake cell, upsampled ---------------
    lake30 = np.load(f"{DATA}/lake_mask.npy")[SUB[0]:SUB[1], SUB[2]:SUB[3]]
    levels = {}
    for s in SCEN:
        me = np.load(f"{OUT}/{s}_maxeta.npy")
        wl = np.where(lake30 & np.isfinite(me), me, np.nan)
        idx = ndimage.distance_transform_edt(~np.isfinite(wl), return_distances=False,
                                             return_indices=True)
        filled = wl[tuple(idx)]
        sl = filled[r0 - SUB[0]:r1 - SUB[0], c0 - SUB[2]:c1 - SUB[2]].astype("float32")
        levels[s] = upto(sl, dem.shape, k)

    # --- lake at 1 m -----------------------------------------------------------
    # LiDAR fills open water with a constant surface value, which gives a far crisper
    # shoreline than upsampling the 30 m mask. Fall back to the 30 m mask only where
    # there is no LiDAR at all.
    lk30 = upto(lake30[r0 - SUB[0]:r1 - SUB[0], c0 - SUB[2]:c1 - SUB[2]].astype("float32"),
                dem.shape, k) > 0.5

    # --- datum tie, taken at the water's edge ---------------------------------
    # Comparing LiDAR to Copernicus over land conflates the datum difference with
    # Copernicus's vegetation/building bias (it is a surface model), which would
    # lower the ground spuriously and over-predict flooding. The water's edge is
    # the clean reference: Copernicus flattens the lake to 345.00 m, and the LiDAR
    # water-surface value is the same physical surface in NZVD2016.
    edge = ndimage.binary_dilation(lk30, np.ones((5, 5))) & np.isfinite(dem) & (~lk30)
    lidar_water = float(np.percentile(dem[edge], 5))
    off = lidar_water - LAKE_LEVEL
    print(f"  LiDAR water surface {lidar_water:.2f} m (NZVD2016) vs Copernicus "
          f"{LAKE_LEVEL:.2f} m -> datum tie {off:+.2f} m")

    lk = (np.abs(dem - lidar_water) < 0.05) | (lk30 & ~np.isfinite(dem))
    lk = ndimage.binary_opening(lk, np.ones((3, 3)))
    print(f"  lake at 1 m: {lk.sum()*MPP*MPP/1e4:.1f} ha in window "
          f"(from LiDAR water surface)")

    ground = dem - off                      # into the model's vertical frame
    nodata = ~np.isfinite(ground)
    ground_f = np.where(nodata, 9999.0, ground)

    # The Hawea dam is the barrier between the lake and the river below it. Enforce
    # the same crest the 30 m model uses, so water reaches the river only by
    # overtopping rather than by leaking through a gap in the terrain model.
    cop = np.load(f"{DATA}/dem.npy")[SUB[0]:SUB[1], SUB[2]:SUB[3]]
    zbl = np.load(f"{DATA}/zb_land.npy")[SUB[0]:SUB[1], SUB[2]:SUB[3]]
    bar30 = (zbl > cop + 0.01)[r0 - SUB[0]:r1 - SUB[0], c0 - SUB[2]:c1 - SUB[2]]
    barrier = upto(bar30.astype("float32"), dem.shape, k) > 0.5
    DAM_CREST = 348.0
    ground_f[barrier] = np.maximum(ground_f[barrier], DAM_CREST)
    print(f"  dam crest enforced on {barrier.sum()*MPP*MPP/1e4:.2f} ha at {DAM_CREST} m")

    floods = {}
    for s in SCEN:
        lvl = levels[s]
        wet = (ground_f < lvl) | lk
        lab, n = ndimage.label(wet)
        keep = set(np.unique(lab[lk])) - {0}
        fl = np.isin(lab, list(keep)) & wet & (~lk) & (~nodata)
        floods[s] = fl
        print(f"  {s}: {fl.sum()*MPP*MPP/1e4:6.2f} ha of lake frontage covered, "
              f"peak level {np.nanmax(np.where(lk, lvl, np.nan)):.2f} m "
              f"({np.nanmax(np.where(lk, lvl, np.nan))-LAKE_LEVEL:+.2f} m)")

    np.savez_compressed(f"{OUT}/township_hires.npz",
                        rgb=np.nan_to_num(rgb, nan=255.0).clip(0, 255).astype("uint8"),
                        dem=dem, lake=lk, nodata=nodata, win=np.array(WIN), mpp=MPP,
                        off=off, **{f"flood_{s}": floods[s] for s in SCEN})
    json.dump({s: round(float(floods[s].sum() * MPP * MPP / 1e4), 2) for s in SCEN},
              open(f"{OUT}/township_hires_ha.json", "w"), indent=2)
    print(f"saved {OUT}/township_hires.npz")


# ---------------------------------------------------------------- figure
def plot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patheffects as pe
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.patches import Patch

    z = np.load(f"{OUT}/township_hires.npz")
    rgb, lk, nod = z["rgb"], z["lake"], z["nodata"]
    x0, y0, x1, y1 = z["win"]; mpp = float(z["mpp"])
    ext = (x0, x1, y0, y1)
    BAND = {"S1_moderate": "#0d366b", "S2_large": "#2a78d6", "S3_extreme": "#9ec5f4"}
    INK, INK2 = "#0b0b0b", "#52514e"

    img = np.transpose(rgb, (1, 2, 0)).astype(float) / 255.0
    blank = img.sum(axis=2) < 0.02          # tiles do not cover open water
    img[blank] = np.array([0.82, 0.87, 0.91])
    aspect = (x1 - x0) / (y1 - y0)
    W = 15.5
    fig, ax = plt.subplots(figsize=(W, W / aspect + 1.0), constrained_layout=True)
    ax.imshow(img, extent=ext, origin="upper", interpolation="bilinear")
    # where there is no LiDAR we cannot say - grey it out rather than imply dry
    ax.imshow(np.where(nod & ~lk, 1.0, np.nan), extent=ext, origin="upper",
              cmap=LinearSegmentedColormap.from_list("n", ["#8e8e8e"] * 2), alpha=.55,
              interpolation="nearest")
    for s in reversed(SCEN):
        fl = z[f"flood_{s}"]
        ax.imshow(np.where(fl, 1.0, np.nan), extent=ext, origin="upper",
                  cmap=LinearSegmentedColormap.from_list("b", [BAND[s]] * 2),
                  alpha=.72, interpolation="nearest")


    ha = json.load(open(f"{OUT}/township_hires_ha.json")) if os.path.exists(
        f"{OUT}/township_hires_ha.json") else {}
    ax.set_xticks([]); ax.set_yticks([])
    ax.plot([x1 - 300, x1 - 100], [y0 + 35] * 2, "-", lw=4, color="white",
            solid_capstyle="butt")
    ax.text(x1 - 200, y0 + 48, "200 m", color="white", fontsize=10, ha="center",
            path_effects=[pe.withStroke(linewidth=2.5, foreground="#00000088")])
    handles = [Patch(fc=BAND[s], alpha=.8,
                     label={"S1_moderate": "covered even in S1 (5×10⁵ m³)",
                            "S2_large": "also covered in S2 (5×10⁶ m³)",
                            "S3_extreme": "also covered in S3 (2×10⁷ m³)"}[s])
               for s in SCEN] + [Patch(fc="#8e8e8e", alpha=.6, label="no LiDAR — not assessed")]
    ax.legend(handles=handles, loc="lower left", fontsize=9.5, framealpha=.94,
              edgecolor="#dcdfdb", ncol=2)
    sub = " · ".join(f"{k.split('_')[0]} {v} ha" for k, v in ha.items())
    ax.set_title("Lake Hāwea township foreshore — what the wave covers\n"
                 f"LINZ 0.1 m aerial imagery over LINZ 1 m LiDAR   ·   {sub}",
                 fontsize=12.5, color=INK)
    fig.savefig(f"{OUT}/figures/09_township_hires.png", dpi=150, facecolor="white")
    plt.close(fig)
    print("wrote 09_township_hires.png")


if __name__ == "__main__":
    import sys
    if "--plot-only" not in sys.argv:
        main()
    plot()
