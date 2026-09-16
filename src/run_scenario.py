"""
Run one landslide-tsunami scenario on Lake Hawea end to end.

Usage:  python3 src/run_scenario.py <scenario> [duration_s]
"""
import sys, os, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from swe_model import ShallowWater
import impulse_wave as iw

DATA = "/home/user/Lake-tsunami-scenarios/data"
OUT  = "/home/user/Lake-tsunami-scenarios/outputs"
RES, LAKE_LEVEL = 30.0, 345.0
XMIN, YMAX = 1292000.0, 5096000.0
SUB = (80, 1520, 100, 990)          # r0, r1, c0, c1 of the model sub-window

# Source zone 4: mid-lake east shore, the highest-ranked EIL zone (38.5 deg, 1125 m relief)
IMPACT_XY = (1308080.0, 5069930.0)

SCENARIOS = {
    # name        volume   width  length  thick  drop  alpha  delta
    "S1_moderate": dict(volume=5.0e5, width=250., length=250., thickness=8.,
                        drop=200., alpha=38.5, delta=25.),
    "S2_large":    dict(volume=5.0e6, width=500., length=400., thickness=25.,
                        drop=300., alpha=38.5, delta=25.),
    "S3_extreme":  dict(volume=2.0e7, width=800., length=700., thickness=36.,
                        drop=400., alpha=38.5, delta=22.),
}

GAUGES = {
    "Hawea township": (1302581, 5053435),
    "Hawea dam":      (1302244, 5052629),
    "The Neck":       (1301033, 5066822),
    "Head of lake":   (1313722, 5086314),
}
GAUGE_MIN_DEPTH = 10.0      # snap to water this deep so a gauge cannot dry out


def rc(x, y):
    return int(round((YMAX - y) / RES)), int(round((x - XMIN) / RES))


def main(name, duration=1800.0, snap_every=10.0):
    sc = SCENARIOS[name]
    dem   = np.load(f"{DATA}/dem.npy")
    lake  = np.load(f"{DATA}/lake_mask.npy")
    depth = np.load(f"{DATA}/depth.npy")

    r0, r1, c0, c1 = SUB
    dem, lake, depth = dem[r0:r1, c0:c1], lake[r0:r1, c0:c1], depth[r0:r1, c0:c1]
    ny, nx = dem.shape

    # bed elevation: lake bed under the lake, DEM surface on land
    # zb_land carries the enforced dam crest at the outlet (see build_bathymetry.py)
    zb_land = np.load(f"{DATA}/zb_land.npy")[r0:r1, c0:c1]
    zb = np.where(lake, LAKE_LEVEL - depth, zb_land).astype(np.float64)
    eta = np.where(lake, LAKE_LEVEL, zb).astype(np.float64)

    # --- source ----------------------------------------------------------
    ir, ic = rc(*IMPACT_XY); ir -= r0; ic -= c0
    yy, xx = np.mgrid[0:ny, 0:nx]
    rad = np.hypot((yy - ir) * RES, (xx - ic) * RES)

    # still-water depth in the slide impact zone: mean over the slide footprint offshore
    foot = lake & (rad < max(sc["width"], 600.0))
    h = float(depth[foot].mean())

    vel = iw.slide_velocity(sc["drop"], sc["alpha"], sc["delta"])
    par = iw.impulse_parameters(sc["volume"], sc["width"], sc["thickness"],
                                vel, h, sc["alpha"])
    wav = iw.near_field_wave(par["P"], h)

    a_M, x_M, L_M = wav["a_M"], wav["x_M"], wav["L_M"]
    Lc = min(np.sqrt(4.0 * h ** 3 / (3.0 * a_M)), L_M / 2.0)

    crest  = a_M * np.cosh(np.clip((rad - x_M) / Lc, -20, 20)) ** -2
    trough = 0.25 * a_M * np.cosh(np.clip((rad - x_M + 1.5 * Lc) / Lc, -20, 20)) ** -2
    disp = np.where(lake, crest - trough, 0.0)

    eta = eta + disp
    m = ShallowWater(zb, eta, RES, manning=0.025, hmin=0.05, cfl=0.35)

    # radially outgoing discharge so the source propagates away from the impact
    with np.errstate(invalid="ignore", divide="ignore"):
        qr = np.where(lake, disp * np.sqrt(9.81 * np.maximum(depth, 1.0)), 0.0)
        ux = np.where(rad > 1, (xx - ic) * RES / np.maximum(rad, 1), 0.0)
        uy = np.where(rad > 1, (yy - ir) * RES / np.maximum(rad, 1), 0.0)
    m.M[:] = qr * ux
    m.N[:] = qr * uy

    print(f"[{name}] volume={sc['volume']:.2e} m3  impact depth h={h:.0f} m")
    print(f"  slide velocity Vs = {vel:.1f} m/s   F={par['F']:.2f} S={par['S']:.3f} "
          f"M={par['M']:.2f}  P={par['P']:.3f}")
    print(f"  near field: H_M={wav['H_M']:.1f} m  a_M={a_M:.1f} m  x_M={x_M:.0f} m  "
          f"T_M={wav['T_M']:.0f} s  L_M={L_M/1000:.1f} km")
    print(f"  dt={m.dt:.3f} s, grid {ny}x{nx}, simulating {duration:.0f} s")

    # --- run --------------------------------------------------------------
    maxeta = np.full_like(eta, -1e9)
    arrival = np.full(eta.shape, np.nan, dtype="float32")
    # snap each gauge to the nearest lake cell so it records water level, not dry ground
    deep = lake & (depth >= GAUGE_MIN_DEPTH)
    ly, lx = np.where(deep)
    gid = {}
    for k, v in GAUGES.items():
        gr, gc = rc(*v); gr -= r0; gc -= c0
        j = int(np.argmin((ly - gr) ** 2 + (lx - gc) ** 2))
        gid[k] = (int(ly[j]), int(lx[j]))
        print(f"  gauge {k:15s} snapped {np.hypot(ly[j]-gr, lx[j]-gc)*RES:5.0f} m, "
              f"depth {depth[gid[k]]:.1f} m")
    series = {k: [] for k in gid}
    times, snaps = [], []
    nsnap = 0
    t0 = time.time()
    nsteps = int(duration / m.dt)
    for k in range(nsteps):
        t = m.step()
        wet = m.eta > m.zb + 0.05
        np.maximum(maxeta, np.where(wet, m.eta, -1e9), out=maxeta)
        first = np.isnan(arrival) & wet & (m.eta - LAKE_LEVEL > 0.5)
        arrival[first] = t
        if t >= nsnap * snap_every:
            snaps.append((np.where(lake, m.eta - LAKE_LEVEL, np.nan)
                          .astype("float32"))[::2, ::2].copy())
            times.append(t); nsnap += 1
            for key, (gr, gc) in gid.items():
                series[key].append(float(m.eta[gr, gc] - LAKE_LEVEL))
    print(f"  ran {nsteps} steps in {time.time()-t0:.0f} s")

    maxeta[maxeta < -1e8] = np.nan
    np.save(f"{OUT}/{name}_maxeta.npy", maxeta.astype("float32"))
    np.save(f"{OUT}/{name}_arrival.npy", arrival)
    np.save(f"{OUT}/{name}_snaps.npy", np.array(snaps))
    meta = dict(name=name, scenario=sc, impact_depth_m=round(h, 1),
                slide_velocity_ms=round(vel, 1),
                F=round(par["F"], 3), S=round(par["S"], 4), M=round(par["M"], 3),
                P=round(par["P"], 3), H_M=round(wav["H_M"], 1), a_M=round(a_M, 1),
                x_M=round(x_M, 0), T_M=round(wav["T_M"], 1), L_M=round(L_M, 0),
                dt=round(m.dt, 4), duration=duration, snap_every=snap_every,
                sub=SUB, res=RES, times=times, gauges=series,
                impact_xy=IMPACT_XY)
    json.dump(meta, open(f"{OUT}/{name}_meta.json", "w"), indent=2)

    # inundation summary on land
    land = ~lake
    inund = land & np.isfinite(maxeta) & (maxeta > zb + 0.1)
    print(f"  inundated land area: {inund.sum()*RES*RES/1e6:.2f} km2")
    for key, (gr, gc) in gid.items():
        s = np.array(series[key])
        print(f"  gauge {key:15s} max {s.max():6.2f} m  min {s.min():6.2f} m")
    return meta


if __name__ == "__main__":
    name = sys.argv[1]
    dur = float(sys.argv[2]) if len(sys.argv) > 2 else 1800.0
    os.makedirs(OUT, exist_ok=True)
    main(name, dur)
