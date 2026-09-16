"""
Build the Lake Hawea model grid: topography, lake mask and synthetic bathymetry.

Source topography : Copernicus GLO-30 DEM (1 arc-second, EGM2008 vertical datum),
                    ESA / Airbus, freely redistributable.
Working projection: NZTM2000 (EPSG:2193), metres, so wave celerity and distances
                    can be handled in a metric grid.
"""
import os, numpy as np, rasterio
from rasterio.merge import merge
from rasterio.warp import calculate_default_transform, reproject, Resampling

SP   = os.environ.get("SCRATCH", "/tmp/sp")
DEM  = os.environ.get("DEM_DIR", f"{SP}/dem")
OUT  = "/home/user/Lake-tsunami-scenarios/data"
os.makedirs(OUT, exist_ok=True)

# Model domain in NZTM2000 (covers Lake Hawea + Hawea township + Albert Town / Wanaka outlet)
NZTM_BOUNDS = (1292000, 5040000, 1324000, 5096000)   # xmin, ymin, xmax, ymax
RES = 30.0                                            # grid resolution, metres

def main():
    tiles = [rasterio.open(f"{DEM}/{f}") for f in sorted(os.listdir(DEM)) if f.endswith(".tif")]
    mosaic, transform = merge(tiles, bounds=(169.0, -44.90, 169.70, -44.20))
    meta = tiles[0].meta.copy()
    meta.update(height=mosaic.shape[1], width=mosaic.shape[2], transform=transform, count=1)
    print(f"mosaic (WGS84): {mosaic.shape}")

    xmin, ymin, xmax, ymax = NZTM_BOUNDS
    width  = int((xmax - xmin) / RES)
    height = int((ymax - ymin) / RES)
    dst_transform = rasterio.transform.from_origin(xmin, ymax, RES, RES)
    dst = np.full((height, width), np.nan, dtype="float32")

    reproject(source=mosaic[0], destination=dst,
              src_transform=transform, src_crs="EPSG:4326",
              dst_transform=dst_transform, dst_crs="EPSG:2193",
              resampling=Resampling.bilinear, src_nodata=None, dst_nodata=np.nan)
    print(f"NZTM grid: {dst.shape} @ {RES} m   z range {np.nanmin(dst):.1f} .. {np.nanmax(dst):.1f} m")

    prof = dict(driver="GTiff", height=height, width=width, count=1, dtype="float32",
                crs="EPSG:2193", transform=dst_transform, nodata=np.nan, compress="deflate")
    with rasterio.open(f"{OUT}/hawea_dem_nztm30.tif", "w", **prof) as f:
        f.write(dst, 1)
    np.save(f"{OUT}/dem.npy", dst)
    with open(f"{OUT}/grid.txt", "w") as f:
        f.write(f"{xmin} {ymin} {xmax} {ymax} {RES} {width} {height}\n")
    return dst

if __name__ == "__main__":
    main()
