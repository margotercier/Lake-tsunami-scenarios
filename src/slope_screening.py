"""
Screen the Lake Hawea shoreline for earthquake-induced landslide (EIL) source areas.

Mirrors the screening logic GNS Science used for lakes Manapouri and Te Anau
(Hancox 2012, GNS CR 2012/146): steep (>=30-45 deg), high-relief slopes that drop
directly into deep water are the credible tsunami-generating sites. Here the screening
is done quantitatively from the DEM instead of by air-photo interpretation.
"""
import numpy as np, json
from scipy import ndimage

DATA = "/home/user/Lake-tsunami-scenarios/data"
RES  = 30.0
LAKE_LEVEL = 345.0
SLOPE_MIN  = 30.0     # deg, threshold for high landslide susceptibility
RELIEF_MIN = 400.0    # m, above lake level within 1.5 km inland


def slope_degrees(dem, res):
    gy, gx = np.gradient(dem, res, res)
    return np.degrees(np.arctan(np.hypot(gx, gy)))


def main():
    dem   = np.load(f"{DATA}/dem.npy")
    lake  = np.load(f"{DATA}/lake_mask.npy")
    depth = np.load(f"{DATA}/depth.npy")
    land  = ~lake

    slope = slope_degrees(dem, RES)

    # inland statistics, computed over land cells only
    lm = land.astype(float)
    r_relief = int(1500 / RES)
    dem_land = np.where(land, dem, -9999.0)
    relief = ndimage.maximum_filter(dem_land, size=2 * r_relief + 1) - LAKE_LEVEL

    r_slope = int(750 / RES)
    num = ndimage.uniform_filter(np.where(land, slope, 0.0), size=2 * r_slope + 1)
    den = ndimage.uniform_filter(lm, size=2 * r_slope + 1)
    mean_slope = np.where(den > 0.05, num / np.maximum(den, 1e-6), 0.0)

    # shoreline = lake cells adjacent to land
    shore = lake & ~ndimage.binary_erosion(lake, np.ones((3, 3)))
    # near-shore water depth 300 m offshore (controls wave generation efficiency)
    d300 = ndimage.maximum_filter(depth, size=2 * int(300 / RES) + 1)

    ys, xs = np.where(shore)
    sc_slope  = mean_slope[ys, xs]
    sc_relief = relief[ys, xs]
    sc_depth  = d300[ys, xs]

    cand = (sc_slope >= SLOPE_MIN) & (sc_relief >= RELIEF_MIN)
    print(f"shoreline cells: {len(ys)}  ({len(ys)*RES/1000:.0f} km of shoreline)")
    print(f"candidate EIL shoreline cells: {cand.sum()} "
          f"({100*cand.sum()/len(ys):.0f}% of shoreline)")
    print(f"  mean slope over candidates : {sc_slope[cand].mean():.1f} deg")
    print(f"  mean relief over candidates: {sc_relief[cand].mean():.0f} m")
    print(f"  mean near-shore depth      : {sc_depth[cand].mean():.0f} m")

    # cluster candidate shoreline cells into contiguous source zones
    cmask = np.zeros_like(lake)
    cmask[ys[cand], xs[cand]] = True
    lab, n = ndimage.label(ndimage.binary_dilation(cmask, np.ones((5, 5))))
    zones = []
    for i in range(1, n + 1):
        sel = (lab == i) & cmask
        k = int(sel.sum())
        if k < 20:                     # ignore segments shorter than ~600 m
            continue
        zy, zx = np.where(sel)
        zones.append(dict(
            id=len(zones) + 1,
            n_cells=k,
            shore_length_km=round(k * RES / 1000, 2),
            row=int(zy.mean()), col=int(zx.mean()),
            mean_slope_deg=round(float(mean_slope[zy, zx].mean()), 1),
            max_relief_m=round(float(relief[zy, zx].max()), 0),
            mean_relief_m=round(float(relief[zy, zx].mean()), 0),
            nearshore_depth_m=round(float(d300[zy, zx].mean()), 0),
        ))
    zones.sort(key=lambda z: z["shore_length_km"] * z["mean_relief_m"], reverse=True)
    print(f"\n{len(zones)} candidate source zones (ranked by length x relief):")
    print(f"{'id':>3} {'len_km':>7} {'slope':>6} {'relief':>7} {'depth':>6}  {'row':>5} {'col':>5}")
    for z in zones[:12]:
        print(f"{z['id']:>3} {z['shore_length_km']:>7.2f} {z['mean_slope_deg']:>6.1f} "
              f"{z['mean_relief_m']:>7.0f} {z['nearshore_depth_m']:>6.0f} "
              f"{z['row']:>5} {z['col']:>5}")

    np.save(f"{DATA}/slope.npy", slope.astype("float32"))
    np.save(f"{DATA}/shore.npy", shore)
    json.dump(zones, open(f"{DATA}/eil_zones.json", "w"), indent=2)

if __name__ == "__main__":
    main()
