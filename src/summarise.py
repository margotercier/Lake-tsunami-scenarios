"""Collect model results into a compact table (used to build RESULTS.md)."""
import json, numpy as np, os
DATA, OUT = "/home/user/Lake-tsunami-scenarios/data", "/home/user/Lake-tsunami-scenarios/outputs"
RES, LAKE_LEVEL, SUB = 30.0, 345.0, (80, 1520, 100, 990)
SCEN = ["S1_moderate", "S2_large", "S3_extreme"]

def sub(a): return a[SUB[0]:SUB[1], SUB[2]:SUB[3]]

dem, lake, depth = (sub(np.load(f"{DATA}/{n}.npy")) for n in ("dem", "lake_mask", "depth"))
zb = np.where(lake, LAKE_LEVEL - depth, dem)
rows = []
for s in SCEN:
    m = json.load(open(f"{OUT}/{s}_meta.json"))
    me = np.load(f"{OUT}/{s}_maxeta.npy")
    wave = np.where(lake & np.isfinite(me), me - LAKE_LEVEL, np.nan)
    inund = (~lake) & np.isfinite(me) & (me > zb + 0.1)
    d = np.where(inund, me - zb, np.nan)
    g = m["gauges"]
    tt = np.array(m["times"])
    town = np.array(g["Hawea township"])
    i = int(np.argmax(np.abs(town) > 0.5)) if (np.abs(town) > 0.5).any() else -1
    rows.append(dict(
        scenario=s, volume=m["scenario"]["volume"], P=m["P"], H_M=m["H_M"], a_M=m["a_M"],
        T_M=m["T_M"], x_M=m["x_M"], h=m["impact_depth_m"], vel=m["slide_velocity_ms"],
        max_wave=float(np.nanmax(wave)),
        p99_wave=float(np.nanpercentile(wave, 99)),
        shore_median=float(np.nanmedian(wave[lake])),
        inund_km2=float(inund.sum() * RES * RES / 1e6),
        inund_max_depth=float(np.nanmax(d)) if inund.any() else 0.0,
        town_max=float(town.max()), town_min=float(town.min()),
        town_arrival_min=float(tt[i] / 60) if i > 0 else None,
        dam_max=float(np.array(g["Hawea dam"]).max()),
    ))

print(f"{'scenario':>12} {'V (m3)':>9} {'P':>5} {'H_M':>6} {'a_M':>6} {'T_M':>5} "
      f"{'maxwave':>8} {'p99':>6} {'inund km2':>10} {'maxdepth':>9} "
      f"{'town+':>6} {'town-':>6} {'arr min':>8} {'dam+':>6}")
for r in rows:
    print(f"{r['scenario']:>12} {r['volume']:9.1e} {r['P']:5.2f} {r['H_M']:6.1f} "
          f"{r['a_M']:6.1f} {r['T_M']:5.0f} {r['max_wave']:8.1f} {r['p99_wave']:6.1f} "
          f"{r['inund_km2']:10.2f} {r['inund_max_depth']:9.1f} {r['town_max']:6.2f} "
          f"{r['town_min']:6.2f} "
          f"{(r['town_arrival_min'] if r['town_arrival_min'] else float('nan')):8.1f} "
          f"{r['dam_max']:6.2f}")
json.dump(rows, open(f"{OUT}/summary.json", "w"), indent=2)

print("\nLake parameters:", json.load(open(f"{DATA}/lake_params.json")))
sm = json.load(open(f"{DATA}/seiche_modes.json"))
print(f"Seiche: Merian T1 = {sm['merian_T1_s']/60:.1f} min; "
      f"2D modes (min): {[round(m['period_min'],1) for m in sm['modes']]}")
