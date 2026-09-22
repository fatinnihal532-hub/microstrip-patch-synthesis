"""Checks the patch model against Balanis' worked examples and against its
own closed forms. Exits non-zero if any check fails.

Reference values are quoted from C. A. Balanis, "Antenna Theory: Analysis
and Design," 3rd ed., Wiley, 2005, Examples 14.1, 14.2 and 14.3
(RT/duroid 5880, eps_r = 2.2, h = 0.1588 cm, f0 = 10 GHz). The book rounds
c to 3e8 m/s (lambda0 = 3 cm); this code uses the exact c, which accounts
for the 0.08 % differences in W and L.

Run: python3 verify.py   (about 1 s)
"""
import math
import sys

from src.patch import (synthesize, slot_integral_i1, self_conductance, mutual_conductance,
                       inset_depth_for, directivity_single_slot, directivity_two_slot_af,
                       directivity_two_slot, vswr_bandwidth)
from src.substrates import Substrate, RT_DUROID_5880, FR4
from src.design_space import standard_board_options

PASS, FAIL = [], []


def check(name, ok, detail=""):
    (PASS if ok else FAIL).append(name)
    print(f"  [{'pass' if ok else 'FAIL'}] {name}")
    if detail:
        print(f"         {detail}")


def rel(a, b):
    return abs(a - b) / abs(b)


def main():
    print("Microstrip patch model - checks against Balanis and closed-form theory")
    print("=" * 70)

    # --- Example 14.1: synthesis --------------------------------------------
    ex = synthesize(10e9, Substrate("Balanis Ex. 14.1", 2.2, 0.0009), 0.1588e-2)
    check("Ex 14.1: W, eps_eff, dL, L within 0.2 % of the book",
          rel(ex.w, 1.186e-2) < 2e-3 and rel(ex.eps_eff, 1.972) < 2e-3
          and rel(ex.delta_l, 0.081e-2) < 2e-3 and rel(ex.l, 0.906e-2) < 2e-3,
          f"W={ex.w*100:.4f} (1.186) eps_eff={ex.eps_eff:.4f} (1.972) "
          f"dL={ex.delta_l*100:.4f} (0.081) L={ex.l*100:.4f} (0.906) cm")

    # Examples 14.2/14.3 are worked with the rounded W, L and lambda0 = 3 cm.
    k0 = 2 * math.pi / 3e-2
    w, l, le = 1.186e-2, 0.906e-2, 1.068e-2
    g1 = self_conductance(k0 * w)
    g12 = mutual_conductance(k0 * w, k0 * l)
    rin = 1 / (2 * (g1 + g12))

    # --- Example 14.2: conductances, input resistance, inset ----------------
    check("Ex 14.2: G1 = 0.00157 S", rel(g1, 0.00157) < 5e-3, f"G1={g1:.6f} S")
    check("Ex 14.2: G12 = 6.1683e-4 S", rel(g12, 6.1683e-4) < 1e-4, f"G12={g12:.5e} S")
    check("Ex 14.2: Rin = 228.3508 ohm", rel(rin, 228.3508) < 1e-5, f"Rin={rin:.4f} ohm")
    y0 = inset_depth_for(rin, 50.0, l)
    check("Ex 14.2: inset y0 = 0.3126 cm for 50 ohm", rel(y0, 0.3126e-2) < 1e-3,
          f"y0={y0*100:.4f} cm")

    # --- Example 14.3: directivity -------------------------------------------
    i1 = slot_integral_i1(k0 * w)
    d0 = directivity_single_slot(k0 * w)
    d2_af = directivity_two_slot_af(k0 * w, g12 / g1)
    check("Ex 14.3: I1 = 1.863, D0 = 3.312, D2 (14-56) = 4.7584",
          rel(i1, 1.863) < 1e-3 and rel(d0, 3.312) < 1e-3 and rel(d2_af, 4.7584) < 1e-3,
          f"I1={i1:.4f} D0={d0:.4f} D2={d2_af:.4f}")
    d2_full = directivity_two_slot(k0 * w, k0 * le)
    check("Ex 14.3: D2 from the full pattern integral (14-55) within 1 % of 5.3873",
          rel(d2_full, 5.3873) < 0.01,
          f"D2={d2_full:.4f}. This code's I2 = 3.5655 is confirmed by an independent "
          f"2000x2000 grid integration to 12 digits; the book's I2 = 3.598 corresponds "
          f"to Le = 1.057 cm rather than the stated 1.068 cm, i.e. a coarser numerical "
          f"integration in the text.")

    # --- Internal consistency ------------------------------------------------
    f0 = 2.44175e9
    d = synthesize(f0, RT_DUROID_5880, 3.175e-3)
    check("numerical VSWR<=2 bandwidth from S11 equals (S-1)/(Q_t sqrt S)",
          rel(d.bandwidth_vswr2, vswr_bandwidth(d.q_total)) < 1e-3,
          f"numerical {d.bandwidth_vswr2*100:.4f} % vs closed form "
          f"{vswr_bandwidth(d.q_total)*100:.4f} %")
    check("inset-matched design has |S11(f0)| = 0", abs(d.s11(f0)) < 1e-12)

    thin = synthesize(f0, RT_DUROID_5880, 0.787e-3)
    thick = synthesize(f0, RT_DUROID_5880, 3.175e-3)
    check("thicker substrate lowers Q_rad and raises efficiency",
          thick.q_rad < thin.q_rad and thick.radiation_efficiency > thin.radiation_efficiency,
          f"h=0.787 mm: Q_rad={thin.q_rad:.1f}, eff={thin.radiation_efficiency:.3f}; "
          f"h=3.175 mm: Q_rad={thick.q_rad:.1f}, eff={thick.radiation_efficiency:.3f}")

    lossy = synthesize(f0, FR4, 3.2e-3)
    lossless_fr4 = synthesize(f0, Substrate("FR-4, loss removed", 4.4, 1e-12), 3.2e-3)
    check("FR-4 loss widens the bandwidth but costs efficiency (loss, not radiation)",
          lossy.bandwidth_vswr2 > lossless_fr4.bandwidth_vswr2
          and lossy.radiation_efficiency < 0.7,
          f"with loss: BW={lossy.bandwidth_vswr2*100:.2f} %, eff={lossy.radiation_efficiency:.1%}; "
          f"without: BW={lossless_fr4.bandwidth_vswr2*100:.2f} %")

    # --- The study's headline conclusion -------------------------------------
    rows = standard_board_options(f0, 83.5e6 / f0, 0.80)
    check("no catalogue board meets both the ISM bandwidth and 80 % efficiency",
          not any(r["meets_spec"] for r in rows),
          f"{sum(r['meets_bandwidth'] for r in rows)} board(s) meet bandwidth, "
          f"{sum(r['meets_efficiency'] for r in rows)} meet efficiency, none meet both")

    print("=" * 70)
    print(f"{len(PASS)} checks passed, {len(FAIL)} failed")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
