"""
Natural seiche (standing-wave) modes of Lake Hawea.

Solves the linear shallow-water eigenvalue problem for a closed basin

        div( h grad phi ) + (omega^2 / g) phi = 0 ,   dphi/dn = 0 on the shore

on the real lake outline and bathymetry. This gives the periods at which the lake
will ring after strong shaking - the quantity that matters for a seismically excited
seiche, and the one a long-period Alpine Fault surface wave train can resonate with.
Compared against the classical Merian approximation T1 = 2L / sqrt(g h).
"""
import numpy as np, json
from scipy import ndimage, sparse
from scipy.sparse.linalg import eigsh

DATA = "/home/user/Lake-tsunami-scenarios/data"
RES, COARSEN = 30.0, 5          # 150 m working grid: seiche wavelengths are kilometres
G = 9.81
NMODES = 8


def main():
    lake = np.load(f"{DATA}/lake_mask.npy")
    depth = np.load(f"{DATA}/depth.npy")

    # coarsen by block-averaging
    ny, nx = lake.shape
    ny2, nx2 = ny // COARSEN, nx // COARSEN
    lk = lake[:ny2 * COARSEN, :nx2 * COARSEN].reshape(ny2, COARSEN, nx2, COARSEN).mean((1, 3))
    dp = depth[:ny2 * COARSEN, :nx2 * COARSEN].reshape(ny2, COARSEN, nx2, COARSEN).mean((1, 3))
    m = lk > 0.5
    dx = RES * COARSEN
    h = np.where(m, np.maximum(dp, 5.0), 0.0)
    print(f"seiche grid {ny2}x{nx2} @ {dx:.0f} m, {m.sum()} wet cells, "
          f"mean depth {h[m].mean():.0f} m")

    idx = -np.ones(m.shape, int)
    idx[m] = np.arange(m.sum())
    N = int(m.sum())

    rows, cols, vals = [], [], []
    diag = np.zeros(N)
    for dy, dx_ in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        nb = np.roll(np.roll(idx, dy, axis=0), dx_, axis=1)
        hn = np.roll(np.roll(h, dy, axis=0), dx_, axis=1)
        both = m & (nb >= 0)
        i = idx[both]; j = nb[both]
        hf = 0.5 * (h[both] + hn[both])          # depth at the shared face
        rows.append(i); cols.append(j); vals.append(hf / dx ** 2)
        np.add.at(diag, i, -hf / dx ** 2)        # no-flux: dry faces contribute nothing
    rows = np.concatenate(rows + [np.arange(N)])
    cols = np.concatenate(cols + [np.arange(N)])
    vals = np.concatenate(vals + [diag])
    L = sparse.csr_matrix((vals, (rows, cols)), shape=(N, N))

    print("solving eigenproblem (shift-invert about 0)...")
    w, v = eigsh(L, k=NMODES + 1, sigma=1e-8, which="LM")
    order = np.argsort(-w)                        # eigenvalues are <= 0; want smallest |w|
    w, v = w[order], v[:, order]

    modes = []
    for k in range(len(w)):
        mu = w[k]
        if mu > -1e-12:                           # the zero mode (uniform level)
            continue
        omega = np.sqrt(-G * mu)
        T = 2 * np.pi / omega
        shape = np.full(m.shape, np.nan)
        shape[m] = v[:, k] / np.abs(v[:, k]).max()
        modes.append(dict(mode=len(modes) + 1, period_s=float(T), period_min=float(T / 60)))
        np.save(f"{DATA}/seiche_mode{len(modes)}.npy", shape.astype("float32"))
        if len(modes) >= NMODES:
            break

    # Merian comparison
    ys, xs = np.where(m)
    Lax = (ys.max() - ys.min()) * dx
    hbar = h[m].mean()
    T_merian = 2 * Lax / np.sqrt(G * hbar)
    print(f"\nlake length (N-S) {Lax/1000:.1f} km, mean depth {hbar:.0f} m")
    print(f"Merian fundamental T1 = 2L/sqrt(gh) = {T_merian:.0f} s = {T_merian/60:.1f} min")
    print(f"\n{'mode':>5} {'period (s)':>11} {'period (min)':>13}")
    for md in modes:
        print(f"{md['mode']:>5} {md['period_s']:>11.0f} {md['period_min']:>13.2f}")

    json.dump(dict(merian_T1_s=round(T_merian, 0), lake_length_km=round(Lax / 1000, 1),
                   mean_depth_m=round(float(hbar), 0), modes=modes),
              open(f"{DATA}/seiche_modes.json", "w"), indent=2)


if __name__ == "__main__":
    main()
