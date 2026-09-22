"""
Substrate / thickness design-space study for a single rectangular patch,
and a design-to-specification solver.

The engineering question this module answers: for a given band and a
given minimum radiation efficiency, which laminate and thickness give the
*smallest* patch that still covers the band at VSWR <= 2 -- and what does
each alternative cost in efficiency, gain and size?

All numbers come from src/patch.py (analytical cavity / transmission-line
model); nothing here introduces new physics.
"""
import math

from scipy.optimize import brentq

from .patch import synthesize, C0
from .substrates import CATALOGUE

# Balanis Sec. 14.1: patches are "usually" designed on 0.003 lambda0 <= h <= 0.05 lambda0.
# The transmission-line/cavity model used here is a thin-substrate model;
# designs outside this window are flagged, not silently trusted.
H_MIN_LAMBDA0 = 0.003
H_MAX_LAMBDA0 = 0.05


def design_row(design):
    return {
        "substrate": design.substrate.name,
        "eps_r": design.substrate.eps_r,
        "tan_delta": design.substrate.tan_delta,
        "h_mm": round(design.h * 1e3, 6),
        "h_over_lambda0": design.h / design.lam0,
        "W_mm": design.w * 1e3,
        "L_mm": design.l * 1e3,
        "area_mm2": design.area_m2 * 1e6,
        "eps_eff": design.eps_eff,
        "R_edge_ohm": design.r_edge,
        "matchable_by_inset": design.matchable_by_inset,
        "inset_mm": design.inset_depth * 1e3,
        "Q_rad": design.q_rad,
        "Q_c": design.q_c,
        "Q_d": design.q_d,
        "Q_t": design.q_total,
        "efficiency": design.radiation_efficiency,
        "bandwidth_pct": design.bandwidth_vswr2 * 100,
        "directivity_dBi": 10 * math.log10(design.directivity),
        "gain_dBi": 10 * math.log10(design.gain),
        "in_thin_substrate_range": H_MIN_LAMBDA0 <= design.h / design.lam0 <= H_MAX_LAMBDA0,
    }


def sweep(f0, thicknesses_m, substrates=CATALOGUE):
    rows = []
    for sub in substrates:
        for h in thicknesses_m:
            rows.append(design_row(synthesize(f0, sub, h)))
    return rows


def minimum_thickness_for_bandwidth(f0, substrate, bw_required, h_lo, h_hi):
    """Smallest h at which the VSWR<=2 fractional bandwidth reaches
    bw_required. Once h is thick enough for radiation rather than copper
    loss to set Q_t, bandwidth rises monotonically with h (Q_rad falls
    roughly as 1/h); every crossing in this study lies in that region, so a
    bracketing root-find is well posed there. Very thin boards can show a
    flat or slightly non-monotonic, loss-dominated bandwidth, which is why
    the bracket starts at 0.003 lambda0. Returns None if even h_hi is not
    enough."""
    def excess(h):
        return synthesize(f0, substrate, h).bandwidth_vswr2 - bw_required
    if excess(h_hi) < 0:
        return None
    if excess(h_lo) >= 0:
        return h_lo
    return brentq(excess, h_lo, h_hi, xtol=1e-7)


def solve_spec(f0, bw_required, min_efficiency, substrates=CATALOGUE):
    """For each laminate: the thinnest board meeting the bandwidth spec
    (inside the thin-substrate validity window), then whether it also
    meets the efficiency spec. Returns one row per laminate with a
    'meets_spec' flag; the recommended design is the smallest-area row
    that meets it."""
    lam0 = C0 / f0
    h_lo, h_hi = H_MIN_LAMBDA0 * lam0, H_MAX_LAMBDA0 * lam0
    out = []
    for sub in substrates:
        h = minimum_thickness_for_bandwidth(f0, sub, bw_required, h_lo, h_hi)
        if h is None:
            out.append({"substrate": sub.name, "meets_spec": False,
                        "reason": "bandwidth not reachable within 0.05 lambda0"})
            continue
        row = design_row(synthesize(f0, sub, h))
        row["meets_spec"] = row["efficiency"] >= min_efficiency
        row["reason"] = "" if row["meets_spec"] else (
            f"efficiency {row['efficiency']:.1%} < {min_efficiency:.0%}: "
            "bandwidth comes from dissipation, not radiation")
        out.append(row)
    return out


def pareto_front(rows, maximize_key, minimize_key):
    """Non-dominated set: maximize one objective, minimize the other."""
    front = []
    for a in rows:
        dominated = any(
            b is not a
            and b[maximize_key] >= a[maximize_key] and b[minimize_key] <= a[minimize_key]
            and (b[maximize_key] > a[maximize_key] or b[minimize_key] < a[minimize_key])
            for b in rows)
        if not dominated:
            front.append(a)
    return sorted(front, key=lambda r: r[minimize_key])


def standard_board_options(f0, bw_required, min_efficiency, substrates=CATALOGUE):
    """Every catalogue board, evaluated -- the realistic version of
    solve_spec(): a designer cannot order a 3.23 mm laminate, only the
    thicknesses a supplier stocks. Each row records whether that board
    meets the bandwidth spec, the efficiency spec, and both."""
    out = []
    for sub in substrates:
        for h_mm in sub.standard_thicknesses_mm:
            row = design_row(synthesize(f0, sub, h_mm * 1e-3))
            row["meets_bandwidth"] = row["bandwidth_pct"] / 100 >= bw_required
            row["meets_efficiency"] = row["efficiency"] >= min_efficiency
            row["meets_spec"] = row["meets_bandwidth"] and row["meets_efficiency"]
            out.append(row)
    return out
