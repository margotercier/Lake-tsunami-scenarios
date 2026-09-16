"""
Landslide-generated impulse wave equations.

Heller V., Hager W.H., Minor H.-E. (2009) "Landslide generated impulse waves in
reservoirs - Basics and computation", VAW-Mitteilung 211, ETH Zurich; 2nd edition
(Evers, Heller, Fuchs, Hager & Boes, 2019). Equation numbers below follow the 2nd
edition. Run-up cross-check uses Synolakis (1987) solitary-wave run-up.

All equations are 2D (channel-shape) generation relations. They are used here only to
size the SOURCE; radial spreading and real basin geometry are handled by the 2D
shallow-water model, which is more appropriate than the idealised 3D relations for a
basin of Lake Hawea's intermediate shape (cf. Ruffini et al. 2019).
"""
import numpy as np

G      = 9.81
RHO_W  = 1000.0          # kg/m3
RHO_S  = 1700.0          # kg/m3, BULK density of a rock avalanche deposit
                         # (~2700 kg/m3 rock at ~37% porosity); keeps D = rho_s/rho_w
                         # inside the 0.59-1.72 range of the VAW-211 experiments


def slide_velocity(drop_height, alpha_deg, delta_deg=25.0):
    """Impact velocity of the slide centre of mass (Korner 1976 energy balance).

    drop_height : vertical fall of the centre of mass to the waterline [m]
    alpha_deg   : slope / impact angle [deg]
    delta_deg   : dynamic basal friction angle [deg]; 20-30 deg for rock avalanches
    """
    a = np.radians(alpha_deg)
    d = np.radians(delta_deg)
    eff = max(0.0, 1.0 - np.tan(d) / np.tan(a))
    return float(np.sqrt(2.0 * G * drop_height * eff))


def impulse_parameters(volume, width, thickness, velocity, depth, alpha_deg,
                       rho_s=RHO_S, rho_w=RHO_W):
    """Dimensionless slide parameters and the impulse product parameter P (Eq. 3.12)."""
    F = velocity / np.sqrt(G * depth)                      # slide Froude number
    S = thickness / depth                                  # relative slide thickness
    M = rho_s * volume / (rho_w * width * depth ** 2)       # relative slide mass
    alpha_eff = np.radians((6.0 / 7.0) * alpha_deg)
    P = F * S ** 0.5 * M ** 0.25 * np.cos(alpha_eff) ** 0.5
    return dict(F=F, S=S, M=M, P=P, D=rho_s / rho_w,
                V=volume / (width * depth ** 2), B=width / depth)


# Validity ranges of the 2D generation equations, VAW-211 2nd ed. Table 3-2
LIMITS = dict(F=(0.86, 6.83), S=(0.09, 1.64), M=(0.11, 10.02), D=(0.59, 1.72),
              V=(0.05, 5.94), B=(0.74, 3.33), P=(0.17, 8.13))
BREAKING_EPS = 0.78      # max relative amplitude of a non-breaking solitary wave


def check_validity(par):
    """Flag governing parameters that fall outside the experimental range."""
    out = {}
    for k, (lo, hi) in LIMITS.items():
        if k in par:
            v = par[k]
            out[k] = "ok" if lo <= v <= hi else ("LOW" if v < lo else "HIGH")
    return out


def near_field_wave(P, depth):
    """Maximum wave in the slide impact zone (Eqs. 3.13-3.18)."""
    H_M = (5.0 / 9.0) * P ** 0.8 * depth                   # (3.13) max wave height
    x_M = (11.0 / 2.0) * P ** 0.5 * depth                  # (3.14) distance to H_M
    T_M = 9.0 * P ** 0.5 * np.sqrt(depth / G)              # (3.15) period of H_M
    a_M = (4.0 / 5.0) * H_M                                # (3.16) crest amplitude
    # A solitary wave cannot exceed ~0.78 h without breaking. Where the empirical
    # relation exceeds that, the near field is a breaking bore: cap the amplitude and
    # record the excess, which is dissipated in the impact zone rather than radiated.
    a_break = BREAKING_EPS * depth
    breaking = a_M > a_break
    a_M_eff = min(a_M, a_break)
    c   = np.sqrt(G * (depth + a_M_eff))                   # (3.17) solitary celerity
    L_M = T_M * c                                          # (3.18) wavelength
    return dict(H_M=H_M, a_M=a_M_eff, a_M_raw=a_M, breaking=bool(breaking),
                x_M=x_M, T_M=T_M, c=c, L_M=L_M)


def far_field_2d(P, depth, x):
    """2D decay in the propagation zone (Eqs. 3.19-3.21), valid for x > x_M."""
    X = np.asarray(x, dtype=float) / depth
    H = (3.0 / 4.0) * P ** 0.8 * X ** (-1.0 / 3.0) * depth       # (3.19)
    T = 9.0 * P ** (5.0 / 16.0) * X ** 0.25 * np.sqrt(depth / G)  # (3.20)
    a = (4.0 / 5.0) * H                                           # (3.16)
    return dict(H=H, a=a, T=T, c=np.sqrt(G * (depth + a)), L=T * np.sqrt(G * (depth + a)))


def runup_synolakis(H, depth, beta_deg):
    """Solitary-wave run-up on a plane slope, Synolakis (1987).

    Non-breaking law  R/h = 2.831 (cot b)^(1/2) (H/h)^(5/4), valid only while
    H/h <= 0.818 (cot b)^(-10/9); beyond that the wave breaks during run-up and the
    law over-predicts badly. Returns (R, is_valid) so callers cannot use it blind.
    This is only an order-of-magnitude cross-check - the reported inundation comes
    from the shallow-water model, which dissipates breaking waves.
    """
    cot = 1.0 / np.tan(np.radians(beta_deg))
    R = 2.831 * np.sqrt(cot) * (H / depth) ** 1.25 * depth
    valid = (H / depth) <= 0.818 * cot ** (-10.0 / 9.0)
    return R, bool(valid)


def breaking_check(a, depth, beta_deg):
    """Grilli et al. (1997) slope parameter; So > 0.37 => non-breaking run-up."""
    eps = a / depth
    return 1.521 * np.tan(np.radians(beta_deg)) / np.sqrt(eps)


# --------------------------------------------------------------------------------
def _validate():
    """Reproduce the Lake Askja 2014 worked example, VAW-211 2nd ed. Section 4.2.1.

    Slide parameters after Ruffini et al. (2019): b = 550 m, s = 35.5 m, alpha = 10.4 deg,
    Vs = 30.1 m/s, m_s = 2e10 kg, h = 138 m.
    Manual reports: P = 0.49, H_M = 43.3 m, x_M = 531 m.
    """
    m_s, b, s, alpha, vel, h = 2e10, 550.0, 35.5, 10.4, 30.1, 138.0
    p = impulse_parameters(m_s / RHO_S, b, s, vel, h, alpha, rho_s=RHO_S)
    w = near_field_wave(p["P"], h)
    print(f"  F = {p['F']:.3f}   (manual 0.82)")
    print(f"  S = {p['S']:.3f}   (manual 0.26)")
    print(f"  M = {p['M']:.3f}   (manual 1.91)")
    print(f"  P = {p['P']:.3f}   (manual 0.49)")
    print(f"  H_M = {w['H_M']:.1f} m   (manual 43.3 m)")
    print(f"  x_M = {w['x_M']:.0f} m   (manual 531 m)")
    ok = (abs(p["P"] - 0.49) < 0.01 and abs(w["H_M"] - 43.3) < 0.8 and abs(w["x_M"] - 531) < 8)
    print("  VALIDATION:", "PASS" if ok else "FAIL")
    return ok


if __name__ == "__main__":
    print("Heller et al. (2009/2019) implementation check - Lake Askja 2014:")
    _validate()
