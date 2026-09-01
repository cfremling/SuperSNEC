#!/usr/bin/env python3
"""Exact Sedov-Taylor blast wave solution for spherical symmetry.

Wraps the Kamm & Timmes (2007) Python implementation (sedov_python)
from https://cococubed.com/research_pages/sedov.shtml

Provides a clean interface: evaluate(r_array, t) -> (rho, vel, press, eps).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

# Add the Timmes sedov_python directory to the path
_SEDOV_DIR = Path(__file__).resolve().parent / "sedov_timmes" / "sedov_python"
if str(_SEDOV_DIR) not in sys.path:
    sys.path.insert(0, str(_SEDOV_DIR))

from globalvars import comvars as gv
from sedov_1d import sed_1d


class SedovSolution:
    """Exact Sedov-Taylor blast wave solution in spherical symmetry."""

    def __init__(self, gamma: float = 1.4, E0: float = 1.0, rho0: float = 1.0):
        self.gamma = gamma
        self.E0 = E0
        self.rho0 = rho0

    def shock_radius(self, t: float) -> float:
        """Compute shock radius by evaluating at one point and reading r2."""
        if t <= 0:
            return 0.0
        gv.gamma = self.gamma
        gv.omega = 0.0
        gv.xgeom = 3.0
        gv.its = 12
        xpos = np.array([0.5])
        sed_1d(t, 1, xpos, self.E0,
               self.rho0, 0.0, 0.0, 0.0, 0.0, gv)
        return gv.r2

    def evaluate(self, r_arr, t: float):
        """Compute exact solution at radii *r_arr* and time *t*.

        Returns (rho, vel, press, eps) arrays with the same shape as r_arr.

        The Timmes solver deletes points where it can't find a root
        (r=0 and any numerically problematic radii), so we feed it
        only positive radii and map results back by matching returned
        xpos to our original grid.
        """
        r_arr = np.atleast_1d(np.asarray(r_arr, dtype=float))
        n = r_arr.size
        den_out = np.zeros(n)
        vel_out = np.zeros(n)
        pres_out = np.zeros(n)
        eps_out = np.zeros(n)

        if t <= 0:
            den_out[:] = self.rho0
            return den_out, vel_out, pres_out, eps_out

        # Only feed positive radii to the solver (it deletes r=0)
        pos = r_arr > 0
        r_pos = r_arr[pos]
        if len(r_pos) == 0:
            return den_out, vel_out, pres_out, eps_out

        # Configure and run
        gv.gamma = self.gamma
        gv.omega = 0.0
        gv.xgeom = 3.0
        gv.its = 12

        xpos = r_pos.copy()
        result = sed_1d(t, len(r_pos), xpos, self.E0,
                        self.rho0, 0.0, 0.0, 0.0, 0.0, gv)
        if result is None:
            return den_out, vel_out, pres_out, eps_out

        den_ret, vel_ret, pres_ret = result[0], result[1], result[2]
        xpos_ret = result[7]

        # The solver deletes near-origin points; interpolate returned
        # profiles back onto our original grid
        if len(xpos_ret) > 1:
            idx_pos = np.where(pos)[0]
            den_out[idx_pos] = np.interp(r_pos, xpos_ret, den_ret,
                                          left=np.nan, right=self.rho0)
            vel_out[idx_pos] = np.interp(r_pos, xpos_ret, vel_ret,
                                          left=np.nan, right=0.0)
            pres_out[idx_pos] = np.interp(r_pos, xpos_ret, pres_ret,
                                           left=np.nan, right=0.0)
            mask = den_out > 0
            eps_out[mask] = pres_out[mask] / ((self.gamma - 1.0) * den_out[mask])

        # Points beyond the shock: rho = rho0 (ambient)
        r2 = gv.r2
        den_out[r_arr >= r2] = self.rho0

        return den_out, vel_out, pres_out, eps_out

    def evaluate_on_grid(self, n_points: int, t: float, r_max: float = None):
        """Return (r, rho, vel, press, eps) on a uniform radial grid."""
        if r_max is None:
            r_max = 1.2 * self.shock_radius(t)
        r_arr = np.linspace(0, r_max, n_points)
        rho, vel, press, eps = self.evaluate(r_arr, t)
        return r_arr, rho, vel, press, eps


# ------------------------------------------------------------------
# Quick self-test
# ------------------------------------------------------------------
if __name__ == "__main__":
    sol = SedovSolution(gamma=1.4, E0=1.0, rho0=1.0)
    t_eval = 4.0
    r2 = sol.shock_radius(t_eval)
    print(f"gamma = {sol.gamma}")
    print(f"t = {t_eval}, r_shock = {r2:.10f}")

    # Profile check
    r_arr, rho, vel, press, eps = sol.evaluate_on_grid(100000, t_eval, r_max=2.0)
    rho_c = np.nan_to_num(rho); vel_c = np.nan_to_num(vel)
    eps_c = np.nan_to_num(eps)
    KE = np.trapezoid(0.5 * rho_c * vel_c**2 * 4 * np.pi * r_arr**2, r_arr)
    TE = np.trapezoid(rho_c * eps_c * 4 * np.pi * r_arr**2, r_arr)
    E_total = KE + TE
    print(f"Energy check: KE={KE:.6e}, TE={TE:.6e}, Total={E_total:.6e} (expect {sol.E0})")
    print(f"Relative energy error: {abs(E_total - sol.E0) / sol.E0:.2e}")

    # Check profiles near center
    print(f"\nProfile at r=0.01: rho={rho[np.argmin(np.abs(r_arr-0.01))]:.4e}, "
          f"vel={vel[np.argmin(np.abs(r_arr-0.01))]:.4e}, "
          f"p={press[np.argmin(np.abs(r_arr-0.01))]:.4e}")
    print(f"Profile at r=0.10: rho={rho[np.argmin(np.abs(r_arr-0.10))]:.4e}, "
          f"vel={vel[np.argmin(np.abs(r_arr-0.10))]:.4e}, "
          f"p={press[np.argmin(np.abs(r_arr-0.10))]:.4e}")
