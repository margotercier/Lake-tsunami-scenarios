"""
Sensitivity of the source wave to the assumed impact-zone water depth h.

VAW-211 explicitly recommends this where h cannot be clearly determined - which is the
case for Lake Hawea because no public bathymetric survey exists. Also reports which
governing parameters fall outside the experimental range of the empirical equations.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import impulse_wave as iw
from run_scenario import SCENARIOS

DEPTHS = [30, 57, 100, 150, 200, 300]

print(f"{'scenario':>12} {'h':>5} {'Vs':>6} {'F':>6} {'S':>6} {'M':>7} {'B':>6} "
      f"{'P':>6} {'H_M':>7} {'a_M':>7} {'brk':>4} {'x_M':>6} {'T_M':>6} {'flags'}")
for name, sc in SCENARIOS.items():
    vel = iw.slide_velocity(sc["drop"], sc["alpha"], sc["delta"])
    for h in DEPTHS:
        par = iw.impulse_parameters(sc["volume"], sc["width"], sc["thickness"],
                                    vel, h, sc["alpha"])
        w = iw.near_field_wave(par["P"], h)
        chk = iw.check_validity(par)
        bad = ",".join(f"{k}={v}" for k, v in chk.items() if v != "ok") or "all in range"
        print(f"{name:>12} {h:5d} {vel:6.1f} {par['F']:6.2f} {par['S']:6.3f} "
              f"{par['M']:7.2f} {par['B']:6.2f} {par['P']:6.2f} {w['H_M']:7.1f} "
              f"{w['a_M']:7.1f} {'Y' if w['breaking'] else 'n':>4} {w['x_M']:6.0f} "
              f"{w['T_M']:6.0f}  {bad}")
    print()

# run-up cross-check (Synolakis 1987); "brk" marks cases where the wave breaks during
# run-up so the non-breaking law does not apply and the model result governs.
print("Run-up cross-check, Synolakis (1987), wave height H at the 50 m depth contour:")
print(f"{'H (m)':>7}" + "".join(f"{b:>12}" for b in ["3 deg", "5 deg", "10 deg", "20 deg"]))
for H in (2, 5, 10, 20):
    cells = []
    for b in (3, 5, 10, 20):
        R, ok = iw.runup_synolakis(H, 50.0, b)
        cells.append(f"{R:8.1f}{'' if ok else ' brk':>4}")
    print(f"{H:>7}" + "".join(cells))
