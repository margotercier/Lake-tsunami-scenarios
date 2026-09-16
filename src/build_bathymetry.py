"""
Lake mask + synthetic bathymetry for Lake Hawea.

No public gridded bathymetry exists for Lake Hawea, so depth is reconstructed from a
distance-to-shore transform calibrated against the published hypsometry
(max depth 392 m, mean depth ~161 m, area ~141 km^2, surface 345 m a.s.l.).
This is the same family of method used by GLOBathy for global lake depth reconstruction.
A longitudinal taper is applied at the northern head to represent Hunter River delta infill.
"""
import numpy as np, rasterio, json
from scipy import ndimage
from scipy.optimize import brentq

DATA = "/home/user/Lake-tsunami-scenarios/data"
LAKE_LEVEL   = 345.0     # m, Copernicus-observed surface (operating range ~338-346 m)
MAX_DEPTH    = 392.0     # m, published
MEAN_DEPTH   = 161.0     # m, published
DELTA_TAPER_KM = 5.0     # Hunter delta infill length at the head

# The Hawea River outlet is controlled by the Hawea dam (Contact Energy). At 30 m the
# Copernicus DEM resolves the outlet channel but not the dam structure across it, so the
# raw DEM leaves the lake hydraulically connected to the river ~7 m below lake level -
# the lake drains through the gap from the first timestep. Enforce a crest across the
# outlet. 348.0 m gives ~3 m of freeboard above the 345 m surface and sits just above the
# 346 m top of the operating range; the true crest should be obtained from the dam owner.
DAM_CREST_M = 348.0

def main():
    dem = np.load(f"{DATA}/dem.npy")
    xmin, ymin, xmax, ymax, res, w, h = open(f"{DATA}/grid.txt").read().split()
    res = float(res)

    # --- lake mask: flat water surface, largest connected component -------------
    water = (dem > LAKE_LEVEL - 0.5) & (dem < LAKE_LEVEL + 0.5)
    water = ndimage.binary_closing(water, np.ones((3, 3)))
    lab, n = ndimage.label(water)
    sizes = ndimage.sum(water, lab, range(1, n + 1))
    lake = lab == (np.argmax(sizes) + 1)
    lake = ndimage.binary_fill_holes(lake)
    area_km2 = lake.sum() * res * res / 1e6
    print(f"lake mask: {lake.sum()} cells, area = {area_km2:.1f} km2 (published 141 km2)")

    # --- distance-to-shore transform ------------------------------------------
    dist = ndimage.distance_transform_edt(lake) * res          # metres from shore
    dnorm = dist / dist.max()
    print(f"max distance to shore: {dist.max()/1000:.2f} km")

    # calibrate exponent p so that mean(depth) matches published mean depth
    vals = dnorm[lake]
    f = lambda p: MAX_DEPTH * np.mean(vals ** p) - MEAN_DEPTH
    p = brentq(f, 0.05, 5.0)
    print(f"calibrated shape exponent p = {p:.3f}  (p<1 => U-shaped glacial trough)")

    depth = np.zeros_like(dem, dtype="float32")
    depth[lake] = MAX_DEPTH * vals ** p

    # --- Hunter delta taper at the northern head ------------------------------
    ys, xs = np.where(lake)
    north_row = ys.min()
    row_dist_km = (np.arange(dem.shape[0]) - north_row) * res / 1000.0
    taper = np.clip(row_dist_km / DELTA_TAPER_KM, 0.0, 1.0) ** 0.7
    depth *= taper[:, None]

    print(f"depth: max {depth.max():.0f} m, mean(lake) {depth[lake].mean():.0f} m, "
          f"volume {depth.sum()*res*res/1e9:.1f} km3")

    # --- seal the outlet so the lake does not drain through the DEM gap ----------
    zb_land = dem.copy()
    below = (~lake) & (dem < LAKE_LEVEL)
    ring = ndimage.binary_dilation(lake, np.ones((3, 3))) & below
    n_leak = int(ring.sum())
    if n_leak:
        barrier = ndimage.binary_dilation(ring, np.ones((3, 3))) & (~lake)
        zb_land[barrier] = np.maximum(zb_land[barrier], DAM_CREST_M)
        print(f"outlet: sealed {n_leak} leak cells "
              f"({barrier.sum()} cells raised to {DAM_CREST_M} m)")
    # verify nothing below lake level is still connected
    below2 = (~lake) & (zb_land < LAKE_LEVEL)
    chk = ndimage.binary_dilation(lake, np.ones((3, 3))) & below2
    grew = int(chk.sum())
    print(f"  remaining connected sub-lake-level cells: {grew} "
          f"({'OK' if grew == 0 else 'STILL LEAKING'})")

    np.save(f"{DATA}/lake_mask.npy", lake)
    np.save(f"{DATA}/depth.npy", depth)
    np.save(f"{DATA}/zb_land.npy", zb_land.astype("float32"))

    prof = dict(driver="GTiff", height=int(h), width=int(w), count=1, dtype="float32",
                crs="EPSG:2193",
                transform=rasterio.transform.from_origin(float(xmin), float(ymax), res, res),
                nodata=0, compress="deflate")
    with rasterio.open(f"{DATA}/hawea_depth_nztm30.tif", "w", **prof) as fo:
        fo.write(depth, 1)

    json.dump(dict(lake_level_m=LAKE_LEVEL, max_depth_m=MAX_DEPTH, mean_depth_m=MEAN_DEPTH,
                   area_km2=round(area_km2, 1), shape_exponent=round(p, 3),
                   volume_km3=round(float(depth.sum()*res*res/1e9), 1), res_m=res),
              open(f"{DATA}/lake_params.json", "w"), indent=2)

if __name__ == "__main__":
    main()
