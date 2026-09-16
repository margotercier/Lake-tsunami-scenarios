"""Generate the Results section of README.md from the model output."""
import json, os
import numpy as np

ROOT = "/home/user/Lake-tsunami-scenarios"
DATA, OUT = f"{ROOT}/data", f"{ROOT}/outputs"
RES, LAKE_LEVEL, SUB = 30.0, 345.0, (80, 1520, 100, 990)
NAMES = {"S1_moderate": "**S1** moderate", "S2_large": "**S2** large",
         "S3_extreme": "**S3** extreme"}
VOLS = {"S1_moderate": "5×10⁵ m³", "S2_large": "5×10⁶ m³", "S3_extreme": "2×10⁷ m³"}


def sub(a):
    return a[SUB[0]:SUB[1], SUB[2]:SUB[3]]


def main():
    rows = json.load(open(f"{OUT}/summary.json"))
    lp = json.load(open(f"{DATA}/lake_params.json"))
    sm = json.load(open(f"{DATA}/seiche_modes.json"))
    zones = json.load(open(f"{DATA}/eil_zones.json"))
    ac = json.load(open(f"{OUT}/areas_covered.json"))
    L = []
    A = L.append

    A("---\n\n## 4. How it was modelled\n")
    A("**Source sizing** — Heller, Hager & Minor (2009/2019) VAW-211 impulse-wave "
      "equations, validated against the manual's own Lake Askja worked example.\n")
    A("**Propagation and inundation** — a 2D nonlinear shallow-water model on a 30 m "
      "NZTM2000 grid (staggered C-grid, upwind face depths, nonlinear advection, "
      "Manning friction, wetting and drying), verified against the analytic long-wave "
      "celerity. The empirical relations size the source only; real basin geometry does "
      "the rest, which is the right call for a lake whose shape is neither the idealised "
      "2D channel nor the 3D open basin of the experiments.\n")
    A(f"**Lake reconstruction** — mask and bathymetry from the DEM's flattened water "
      f"surface at {lp['lake_level_m']:.0f} m a.s.l.: {lp['area_km2']} km² "
      f"(published ~141 km²), max depth {lp['max_depth_m']:.0f} m, mean "
      f"{lp['mean_depth_m']:.0f} m, volume {lp['volume_km3']} km³.\n")

    A("### Where the landslides would come from\n")
    A("Screening the whole 137 km shoreline for slopes ≥30° with ≥400 m of relief "
      "within 1.5 km inland leaves **four candidate source zones**, about 7% of the "
      "shoreline:\n")
    A("| Zone | Position | Shore length | Mean slope | Relief | Distance to township |")
    A("|---|---|---|---|---|---|")
    pos = {4: "mid-lake, east shore", 1: "head of lake (Hunter arm)",
           3: "mid-lake, west shore", 2: "north-east shore"}
    for z in zones:
        A(f"| {z['id']} | {pos.get(z['id'], '—')} ({z['lat']:.3f}°S, {z['lon']:.3f}°E) | "
          f"{z['shore_length_km']:.1f} km | {z['mean_slope_deg']:.0f}° | "
          f"{z['mean_relief_m']:.0f} m | {z['dist_town_km']:.0f} km |")
    A("\nScenarios S1–S3 all place the slide in **zone 4**, the highest-ranked "
      "(38.5° mean slope, 1125 m of relief, 17 km up-lake from the township).\n")

    A("## 5. Results\n")
    A("### Landslide scenarios\n")
    A("| | Volume | P | Near-field H_M | Peak on lake | At township | First arrival | "
      "Largest wave at | At dam | Dam overtops | Land flooded |")
    A("|---|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        arr = (f"{r['town_arrival_min']:.0f} min" if r["town_arrival_min"] else "—")
        ot = f"**+{r['dam_overtop_m']:.2f} m**" if r["overtops"] else "no"
        A(f"| {NAMES[r['scenario']]} | {VOLS[r['scenario']]} | {r['P']:.2f} | "
          f"{r['H_M']:.0f} m | {r['max_wave']:.0f} m | **{r['town_max']:.1f} m** | "
          f"{arr} | {r['town_peak_min']:.0f} min | {r['dam_max']:.1f} m | {ot} | "
          f"{r['inund_km2']:.2f} km² |")
    A("")
    A("Two results matter more than the headline heights.\n")
    A("**The first wave is not the biggest.** At the township the leading wave arrives in "
      f"{min(r['town_arrival_min'] for r in rows if r['town_arrival_min']):.0f}–"
      f"{max(r['town_arrival_min'] for r in rows if r['town_arrival_min']):.0f} minutes at "
      "only 1–1.5 m, and the level then oscillates with *growing* amplitude, peaking "
      f"{rows[1]['town_peak_min']:.0f}–{rows[0]['town_peak_min']:.0f} minutes after the "
      "slide. This is energy being redistributed into the southern basin, not numerical "
      "instability: lake-wide wave energy decays monotonically to about a quarter of its "
      "initial value over the same period, and the largest amplitude anywhere is at t = 0. "
      "For evacuation planning this inverts the usual instinct — the shore is *more* "
      "dangerous twenty minutes in than at first arrival.\n")
    dam = [r for r in rows if r["overtops"]]
    if dam:
        A("**The dam overtops in the larger scenarios.** With the assumed 348.0 m crest "
          "(3 m of freeboard), " + " and ".join(NAMES[r["scenario"]] for r in dam) +
          " overtop, by " + " and ".join(f"{r['dam_overtop_m']:.2f} m" for r in dam) +
          ". This is exactly the concern AF8 raises for Lake Hāwea. **It depends entirely "
          "on the assumed crest elevation**, which is not public — the real figure from the "
          "dam owner could move this result either way, and it is the single cheapest input "
          "that would sharpen this study.\n")

    if os.path.exists(f"{OUT}/S4_tectonic_meta.json"):
        m4 = json.load(open(f"{OUT}/S4_tectonic_meta.json"))
        t4 = np.array(m4["gauges"]["Hawea township"])
        d4 = np.array(m4["gauges"]["Hawea dam"])
        me = np.load(f"{OUT}/S4_tectonic_maxeta.npy")
        dem, lake, depth = (sub(np.load(f"{DATA}/{n}.npy")) for n in
                            ("dem", "lake_mask", "depth"))
        wave = np.where(lake & np.isfinite(me), me - LAKE_LEVEL, np.nan)
        A("### Fault scenario\n")
        A(f"**S4 — Hunter Valley Fault**, the source AF8 names for this lake. Mapped "
          f"trace from NZAFD, strike {m4['strike_deg']:+.0f}°, rupture span "
          f"{m4['rupture_km']:.0f} km, lakebed displacement {m4['uz_max']:+.1f} m "
          f"(hanging wall) to {m4['uz_min']:+.1f} m (footwall).\n")
        A(f"- Peak wave on the lake: **{np.nanmax(wave):.1f} m**")
        A(f"- Lake Hāwea township: **{t4.max():+.1f} m / {t4.min():+.1f} m**")
        A(f"- Hāwea dam: **{d4.max():+.1f} m**\n")
        A("The tectonic wave is smaller than the big landslide cases but it is "
          "*lake-wide and immediate* — the whole water surface is displaced at once, "
          "so there is no travel time to speak of at the southern end, and it is "
          "followed directly by seiching.\n")

    A("### Seiche\n")
    A(f"After the first wave the basin rings. Solving the eigenvalue problem on the "
      f"real outline and bathymetry gives a fundamental period of "
      f"**{sm['modes'][0]['period_min']:.0f} minutes** "
      f"(Merian approximation {sm['merian_T1_s']/60:.0f} min), with higher modes at "
      + ", ".join(f"{m['period_min']:.0f}" for m in sm['modes'][1:5]) + " minutes.\n")
    A("This matters for response: **the hazard does not end when the first wave has "
      "passed.** Water keeps sloshing at these periods for hours, and the second or "
      "third oscillation can be as dangerous as the first for anyone who has gone back "
      "to the shore.\n")

    A("### What gets covered\n")
    A("Most of Hāwea's shoreline is too steep for water to travel far inland, so the "
      "flooded strip is narrow almost everywhere. The exceptions are the deltas — and "
      "the deltas are where the flat, usable land is.\n")
    A("- **Lake Hāwea township** sits on a terrace above the lake. What floods is the "
      "foreshore, the boat ramp and the outlet channel below the dam — "
      f"{ac['township']['S1']:.2f} km² in S1 rising to {ac['township']['S3']:.2f} km² "
      "in S3, not the town itself.")
    A("- **The head of the lake (Hunter River delta)** is the largest area covered "
      f"anywhere on Hāwea: shallow water spreading over {ac['head']['S3']:.2f} km² "
      "of flats in S3.")
    A("- **The slide zone itself** takes run-up of tens of metres onto steep ground — "
      "little area, total destruction.\n")
    A("### Maps and animation\n")
    figs = [("08_areas_covered.png",
             "What actually gets covered — nested flood extents at the places that matter"),
            ("01_overview.png", "Maximum wave field, all three landslide scenarios"),
            ("02_township.png", "Lake Hāwea township inundation detail"),
            ("03_arrival.png", "Wave arrival time"),
            ("04_gauges.png", "Water-level time series"),
            ("05_seiche.png", "Natural seiche modes of the basin"),
            ("06_volume_vs_height.png", "Scenarios against the observed record"),
            ("07_tectonic.png", "Hunter Valley Fault rupture scenario")]
    for f, cap in figs:
        if os.path.exists(f"{OUT}/figures/{f}"):
            A(f"**{cap}**\n\n![{cap}](outputs/figures/{f})\n")
    for v in ("S2_large", "S3_extreme"):
        if os.path.exists(f"{OUT}/{v}_animation.mp4"):
            A(f"Animation — `outputs/{v}_animation.mp4`\n")

    txt = "\n".join(L)
    rd = open(f"{ROOT}/README.template.md").read()
    rd = rd.replace("<!--RESULTS-->", txt)
    open(f"{ROOT}/README.md", "w").write(rd)
    print(f"README.md updated ({len(txt)} chars of results)")


if __name__ == "__main__":
    main()
