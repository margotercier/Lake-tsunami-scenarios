"""
Scenario S4: co-seismic tsunami from rupture of the Hunter Valley Fault beneath Lake Hawea.

The AF8 Hazard Scenario (AF8 2016, p.24) names this fault explicitly: "the Hunter Valley
Fault under Lake Hawea, and a tsunami there would likely impact Lake Hawea township and
possibly the hydro control dam". The trace used here is the mapped one from the GNS
Science New Zealand Active Faults Database (1:250,000), queried live.

Source model: a first-order elastic dislocation - uplift of the hanging wall, lesser
subsidence of the footwall, smoothed across strike - rather than a full Okada (1985)
solution, because the fault's dip, width and slip are all listed as Unknown in NZAFD.
Displacement is applied instantaneously to the free surface, the standard tsunami
initialisation for a fast rupture.

Rupture size after Wells & Coppersmith (1994) for a ~30 km reverse rupture:
Mw ~ 6.8, average slip ~0.7 m, maximum ~1.4 m; vertical component ~0.7 x slip at 45 deg dip.
"""
import os, sys, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from swe_model import ShallowWater

DATA = "/home/user/Lake-tsunami-scenarios/data"
OUT  = "/home/user/Lake-tsunami-scenarios/outputs"
RES, LAKE_LEVEL, XMIN, YMAX = 30.0, 345.0, 1292000.0, 5096000.0
SUB = (80, 1520, 100, 990)
NAME = "S4_tectonic"

UPLIFT_M   = 1.0      # max hanging-wall (east) uplift, m
SUBSIDE_M  = 0.3      # max footwall (west) subsidence, m
WIDTH_M    = 1500.0   # across-strike smoothing half-width, m
TAPER_M    = 3000.0   # along-strike taper beyond the mapped rupture, m

GAUGES = {"Hawea township": (1302581, 5053435), "Hawea dam": (1302244, 5052629)}


def rc(x, y):
    return int(round((YMAX - y) / RES)), int(round((x - XMIN) / RES))


def main(duration=1500.0, snap_every=8.0):
    from pyproj import Transformer
    faults = json.load(open(f"{DATA}/faults_near_hawea.json"))
    pts = np.array(faults["Hunter Valley Fault"]["pts"])
    tr = Transformer.from_crs(4326, 2193, always_xy=True)
    fx, fy = tr.transform(pts[:, 0], pts[:, 1])
    fx, fy = np.asarray(fx), np.asarray(fy)

    dem   = np.load(f"{DATA}/dem.npy")
    lake  = np.load(f"{DATA}/lake_mask.npy")
    depth = np.load(f"{DATA}/depth.npy")
    r0, r1, c0, c1 = SUB
    dem, lake, depth = dem[r0:r1, c0:c1], lake[r0:r1, c0:c1], depth[r0:r1, c0:c1]
    ny, nx = dem.shape

    X = XMIN + (np.arange(nx) + c0) * RES
    Y = YMAX - (np.arange(ny) + r0) * RES
    XX, YY = np.meshgrid(X, Y)

    # straight-line fit to the mapped trace: X = a*Y + b  (the trace is near-linear here)
    a, b = np.polyfit(fy, fx, 1)
    strike = np.degrees(np.arctan2(a, 1.0))
    # signed perpendicular distance, positive to the east (hanging wall of a reverse fault)
    dperp = (XX - (a * YY + b)) / np.sqrt(1.0 + a * a)
    # along-strike coordinate, tapered outside the mapped rupture
    s = (YY + a * XX) / np.sqrt(1.0 + a * a)
    s_f = (fy + a * fx) / np.sqrt(1.0 + a * a)
    s0, s1 = s_f.min(), s_f.max()
    taper = np.clip(np.minimum(s - (s0 - TAPER_M), (s1 + TAPER_M) - s) / TAPER_M, 0.0, 1.0)

    th = np.tanh(dperp / WIDTH_M)
    uz = (UPLIFT_M * 0.5 * (1 + th) - SUBSIDE_M * 0.5 * (1 - th)) * taper
    uz -= np.average(uz, weights=lake) * 0.0     # keep absolute displacement, not demeaned

    zb = np.where(lake, LAKE_LEVEL - depth, dem).astype(np.float64)
    # the lake bed moves with the ground; the free surface inherits that displacement
    zb_new = zb + np.where(lake, uz, 0.0)
    eta = np.where(lake, LAKE_LEVEL + uz, zb_new).astype(np.float64)

    print(f"[{NAME}] Hunter Valley Fault, mapped length {np.hypot(np.diff(fx), np.diff(fy)).sum()/1000:.1f} km")
    print(f"  strike {strike:+.0f} deg from N, rupture span {(s1-s0)/1000:.1f} km")
    print(f"  displacement over the lake: max {uz[lake].max():+.2f} m, min {uz[lake].min():+.2f} m")
    print(f"  net volume displaced: {(uz*lake).sum()*RES*RES/1e6:.2f} million m3")

    m = ShallowWater(zb_new, eta, RES, manning=0.025, hmin=0.05, cfl=0.35)
    print(f"  dt={m.dt:.3f} s, grid {ny}x{nx}, simulating {duration:.0f} s")

    maxeta = np.full_like(eta, -1e9)
    arrival = np.full(eta.shape, np.nan, dtype="float32")
    ly, lx = np.where(lake)
    gid = {}
    for k, v in GAUGES.items():
        gr, gc = rc(*v); gr -= r0; gc -= c0
        j = int(np.argmin((ly - gr) ** 2 + (lx - gc) ** 2))
        gid[k] = (int(ly[j]), int(lx[j]))
    series = {k: [] for k in gid}
    times, snaps, nsnap = [], [], 0
    t0 = time.time()
    for k in range(int(duration / m.dt)):
        t = m.step()
        wet = m.eta > m.zb + 0.05
        np.maximum(maxeta, np.where(wet, m.eta, -1e9), out=maxeta)
        first = np.isnan(arrival) & wet & (m.eta - LAKE_LEVEL > 0.2)
        arrival[first] = t
        if t >= nsnap * snap_every:
            snaps.append(np.where(lake, m.eta - LAKE_LEVEL, np.nan).astype("float32")[::2, ::2].copy())
            times.append(t); nsnap += 1
            for key, (gr, gc) in gid.items():
                series[key].append(float(m.eta[gr, gc] - LAKE_LEVEL))
    print(f"  ran in {time.time()-t0:.0f} s")

    maxeta[maxeta < -1e8] = np.nan
    np.save(f"{OUT}/{NAME}_maxeta.npy", maxeta.astype("float32"))
    np.save(f"{OUT}/{NAME}_arrival.npy", arrival)
    np.save(f"{OUT}/{NAME}_snaps.npy", np.array(snaps))
    json.dump(dict(name=NAME, kind="tectonic", fault="Hunter Valley Fault",
                   uplift_m=UPLIFT_M, subsidence_m=SUBSIDE_M, width_m=WIDTH_M,
                   strike_deg=round(float(strike), 1),
                   rupture_km=round(float((s1 - s0) / 1000), 1),
                   uz_max=round(float(uz[lake].max()), 3), uz_min=round(float(uz[lake].min()), 3),
                   dt=round(m.dt, 4), duration=duration, snap_every=snap_every,
                   sub=SUB, res=RES, times=times, gauges=series),
              open(f"{OUT}/{NAME}_meta.json", "w"), indent=2)

    land = ~lake
    inund = land & np.isfinite(maxeta) & (maxeta > zb_new + 0.1)
    print(f"  inundated land area: {inund.sum()*RES*RES/1e6:.2f} km2")
    for key in gid:
        s_ = np.array(series[key])
        print(f"  gauge {key:15s} max {s_.max():6.2f} m  min {s_.min():6.2f} m")


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    main(float(sys.argv[1]) if len(sys.argv) > 1 else 1500.0)
