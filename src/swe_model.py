"""
Nonlinear shallow-water model for impulse-wave propagation and inundation on Lake Hawea.

Explicit staggered (Arakawa C) finite-difference scheme with upwind depth at faces,
nonlinear advection, Manning bed friction and wetting/drying - the classic formulation
used in operational tsunami codes (Kowalik & Murty 1993; Imamura et al. 2006).

State: eta  free-surface elevation at cell centres [m]
       M,N  depth-integrated discharge per unit width at x- and y-faces [m2/s]
"""
import numpy as np

G = 9.81


class ShallowWater:
    def __init__(self, zb, eta0, dx, manning=0.025, hmin=0.05, cfl=0.35):
        self.zb = zb.astype(np.float64)
        self.eta = eta0.astype(np.float64)
        self.dx = float(dx)
        self.n2 = manning ** 2
        self.hmin = hmin
        self.M = np.zeros_like(self.zb)          # x-face discharge, face i-1/2 stored at i
        self.N = np.zeros_like(self.zb)          # y-face discharge
        hmax = float(np.max(self.eta - self.zb))
        self.dt = cfl * self.dx / (np.sqrt(G * max(hmax, 1.0)) * np.sqrt(2.0))
        self.t = 0.0

    # ---- helpers ---------------------------------------------------------
    def depth(self):
        return np.maximum(self.eta - self.zb, 0.0)

    def _face_depth_x(self):
        """Upwind/total depth at x-faces: max(eta) - max(zb) across the face."""
        eta_l, eta_r = self.eta[:, :-1], self.eta[:, 1:]
        zb_l, zb_r = self.zb[:, :-1], self.zb[:, 1:]
        return np.maximum(np.maximum(eta_l, eta_r) - np.maximum(zb_l, zb_r), 0.0)

    def _face_depth_y(self):
        eta_u, eta_d = self.eta[:-1, :], self.eta[1:, :]
        zb_u, zb_d = self.zb[:-1, :], self.zb[1:, :]
        return np.maximum(np.maximum(eta_u, eta_d) - np.maximum(zb_u, zb_d), 0.0)

    # ---- one time step ---------------------------------------------------
    def step(self):
        dx, dt = self.dx, self.dt

        Dx = self._face_depth_x()                       # (ny, nx-1)
        Dy = self._face_depth_y()                       # (ny-1, nx)
        wetx, wety = Dx > self.hmin, Dy > self.hmin

        # --- momentum: pressure gradient -------------------------------
        dedx = (self.eta[:, 1:] - self.eta[:, :-1]) / dx
        dedy = (self.eta[1:, :] - self.eta[:-1, :]) / dx
        Mf = self.M[:, :-1]
        Nf = self.N[:-1, :]
        Mn = np.where(wetx, Mf - dt * G * Dx * dedx, 0.0)
        Nn = np.where(wety, Nf - dt * G * Dy * dedy, 0.0)

        # --- nonlinear advection (upwind, first order) -------------------
        with np.errstate(divide="ignore", invalid="ignore"):
            ux = np.where(wetx, Mf / np.maximum(Dx, self.hmin), 0.0)
            vy = np.where(wety, Nf / np.maximum(Dy, self.hmin), 0.0)
        dMdx = np.zeros_like(Mn); dNdy = np.zeros_like(Nn)
        dMdx[:, 1:-1] = np.where(ux[:, 1:-1] > 0,
                                 (Mf[:, 1:-1] - Mf[:, :-2]), (Mf[:, 2:] - Mf[:, 1:-1])) / dx
        dNdy[1:-1, :] = np.where(vy[1:-1, :] > 0,
                                 (Nf[1:-1, :] - Nf[:-2, :]), (Nf[2:, :] - Nf[1:-1, :])) / dx
        Mn -= dt * ux * dMdx
        Nn -= dt * vy * dNdy

        # --- Manning friction, semi-implicit -----------------------------
        with np.errstate(divide="ignore", invalid="ignore"):
            fx = dt * G * self.n2 * np.abs(Mn) / np.maximum(Dx, self.hmin) ** (7.0 / 3.0)
            fy = dt * G * self.n2 * np.abs(Nn) / np.maximum(Dy, self.hmin) ** (7.0 / 3.0)
        Mn = np.where(wetx, Mn / (1.0 + fx), 0.0)
        Nn = np.where(wety, Nn / (1.0 + fy), 0.0)

        self.M[:, :-1] = Mn
        self.N[:-1, :] = Nn

        # --- continuity ---------------------------------------------------
        dM = np.zeros_like(self.eta); dN = np.zeros_like(self.eta)
        dM[:, 1:-1] = (self.M[:, 1:-1] - self.M[:, :-2]) / dx
        dN[1:-1, :] = (self.N[1:-1, :] - self.N[:-2, :]) / dx
        self.eta -= dt * (dM + dN)

        # keep dry land at bed level (no negative water depth)
        dry = self.eta < self.zb
        self.eta[dry] = self.zb[dry]
        self.t += dt
        return self.t
