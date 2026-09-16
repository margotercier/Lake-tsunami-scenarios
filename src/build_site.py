"""
Build site/index.html from site/template.html by injecting the model output.

Keeps the shareable page and the repository in lockstep: every number on the page
comes from outputs/, never hand-typed.
"""
import json, os, shutil
import numpy as np

ROOT = "/home/user/Lake-tsunami-scenarios"
OUT, DATA, SITE = f"{ROOT}/outputs", f"{ROOT}/data", f"{ROOT}/site"
LAKE_LEVEL, SUB = 345.0, (80, 1520, 100, 990)
LABELS = {"S1_moderate": "S1 moderate", "S2_large": "S2 large", "S3_extreme": "S3 extreme"}
VOLS = {"S1_moderate": "5×10⁵ m³", "S2_large": "5×10⁶ m³", "S3_extreme": "2×10⁷ m³"}


def main():
    rows = json.load(open(f"{OUT}/summary.json"))
    sm = json.load(open(f"{DATA}/seiche_modes.json"))
    by = {r["scenario"]: r for r in rows}

    payload = dict(
        scenarios=[dict(label=LABELS[r["scenario"]], volume=VOLS[r["scenario"]],
                        P=round(r["P"], 2), H_M=round(r["H_M"]),
                        max_wave=round(r["max_wave"], 1),
                        town_max=round(r["town_max"], 1),
                        arrival=(f"{r['town_arrival_min']:.0f} min"
                                 if r["town_arrival_min"] else "—"),
                        dam_max=round(r["dam_max"], 1),
                        peak_min=round(r["town_peak_min"]),
                        overtop=(f"+{r['dam_overtop_m']:.2f} m" if r["overtops"] else "no"),
                        inund=round(r["inund_km2"], 2)) for r in rows],
        town_peak_s2=round(by["S2_large"]["town_max"], 1),
        arrival_s2=(f"{by['S2_large']['town_arrival_min']:.0f}"
                    if by["S2_large"]["town_arrival_min"] else "—"),
        seiche1=round(sm["modes"][0]["period_min"]),
        seiche_rest=", ".join(f"{m['period_min']:.0f}" for m in sm["modes"][1:5]),
        merian=round(sm["merian_T1_s"] / 60),
        arrival_lo=round(min(r["town_arrival_min"] for r in rows if r["town_arrival_min"])),
        arrival_hi=round(max(r["town_arrival_min"] for r in rows if r["town_arrival_min"])),
        peak_lo=round(min(r["town_peak_min"] for r in rows)),
        peak_hi=round(max(r["town_peak_min"] for r in rows)),
        overtoppers=[r["scenario"].split("_")[0] for r in rows if r["overtops"]],
    )

    if os.path.exists(f"{OUT}/township_hires_ha.json"):
        payload["township_ha"] = json.load(open(f"{OUT}/township_hires_ha.json"))

    if os.path.exists(f"{OUT}/S2_large_animation.mp4"):
        payload["video"] = "S2_large_animation.mp4"

    p4 = f"{OUT}/S4_tectonic_meta.json"
    if os.path.exists(p4):
        m4 = json.load(open(p4))
        me = np.load(f"{OUT}/S4_tectonic_maxeta.npy")
        lake = np.load(f"{DATA}/lake_mask.npy")[SUB[0]:SUB[1], SUB[2]:SUB[3]]
        wave = np.where(lake & np.isfinite(me), me - LAKE_LEVEL, np.nan)
        payload["tectonic"] = dict(uz_max=round(m4["uz_max"], 1),
                                   max_wave=round(float(np.nanmax(wave)), 1),
                                   town=round(float(max(m4["gauges"]["Hawea township"])), 1))

    html = open(f"{SITE}/template.html").read()
    inject = ("<script>window.__HAWEA__ = "
              + json.dumps(payload, ensure_ascii=False) + ";</script>\n\n")
    marker = "<script>\n/* Model output is injected here"
    assert marker in html, "injection marker not found in template"
    html = html.replace(marker, inject + marker)
    open(f"{SITE}/index.html", "w").write(html)

    # figures alongside the page
    dst = f"{SITE}/figures"
    os.makedirs(dst, exist_ok=True)
    for f in sorted(os.listdir(f"{OUT}/figures")):
        if f.endswith(".png"):
            shutil.copy(f"{OUT}/figures/{f}", f"{dst}/{f}")
    for v in ("S2_large", "S3_extreme"):
        src = f"{OUT}/{v}_animation.mp4"
        if os.path.exists(src):
            shutil.copy(src, f"{SITE}/{v}_animation.mp4")
    print(f"site/index.html built ({len(html)/1024:.0f} KB), "
          f"{len(os.listdir(dst))} figures copied")


if __name__ == "__main__":
    main()
