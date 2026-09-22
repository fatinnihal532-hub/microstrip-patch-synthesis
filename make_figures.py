"""Regenerates every figure and CSV in results/ from the model in src/.
Run: python3 make_figures.py   (about 10 s)
"""
import csv
import math
import os

import matplotlib.pyplot as plt
import numpy as np

from src.patch import synthesize, e_plane_pattern, h_plane_pattern
from src.substrates import CATALOGUE, RT_DUROID_5880, FR4, RO4003C
from src.design_space import (sweep, solve_spec, standard_board_options,
                              H_MIN_LAMBDA0, H_MAX_LAMBDA0)
from src.plotstyle import save_light_dark

# The 2.4 GHz ISM band (2400-2483.5 MHz), designed at its centre.
BAND_LO, BAND_HI = 2.400e9, 2.4835e9
F0 = 0.5 * (BAND_LO + BAND_HI)
BW_REQUIRED = (BAND_HI - BAND_LO) / F0
MIN_EFFICIENCY = 0.80
LAMBDA0 = 299_792_458.0 / F0
H_SWEEP = np.linspace(H_MIN_LAMBDA0 * LAMBDA0, H_MAX_LAMBDA0 * LAMBDA0, 60)


def _mark_catalogue(ax, theme, i, sub, ykey, rows_std):
    pts = [r for r in rows_std if r["substrate"] == sub.name]
    ax.scatter([r["h_mm"] for r in pts], [r[ykey] for r in pts], s=22,
               color=theme["series"][i], edgecolor=theme["surface"], linewidth=0.8, zorder=3)


def fig_bandwidth_efficiency_vs_h(theme, rows, rows_std):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 3.8))
    for i, sub in enumerate(CATALOGUE):
        rs = [r for r in rows if r["substrate"] == sub.name]
        h = [r["h_mm"] for r in rs]
        ax1.plot(h, [r["bandwidth_pct"] for r in rs], color=theme["series"][i], lw=1.8, label=sub.name)
        ax2.plot(h, [r["efficiency"] * 100 for r in rs], color=theme["series"][i], lw=1.8, label=sub.name)
        _mark_catalogue(ax1, theme, i, sub, "bandwidth_pct", rows_std)
        pts = [dict(r, eff_pct=r["efficiency"] * 100) for r in rows_std]
        _mark_catalogue(ax2, theme, i, sub, "eff_pct", pts)
    ax1.axhline(BW_REQUIRED * 100, color=theme["secondary"], ls=":", lw=1.2)
    ax1.text(0.1, BW_REQUIRED * 100 + 0.08, "full ISM band", color=theme["secondary"], fontsize=8)
    ax2.axhline(MIN_EFFICIENCY * 100, color=theme["secondary"], ls=":", lw=1.2)
    ax2.text(0.1, MIN_EFFICIENCY * 100 + 1.5, "80 % efficiency", color=theme["secondary"], fontsize=8)
    ax1.set_xlabel("Substrate thickness h (mm)")
    ax1.set_ylabel("VSWR<=2 bandwidth (%)")
    ax2.set_xlabel("Substrate thickness h (mm)")
    ax2.set_ylabel("Radiation efficiency (%)")
    ax1.set_title("Bandwidth vs. thickness")
    ax2.set_title("Efficiency vs. thickness")
    ax1.legend(frameon=False, fontsize=7.5, loc="upper left")
    fig.suptitle("Dots: catalogue thicknesses a board can actually be bought in", fontsize=8.5,
                 color=theme["secondary"], y=0.02)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    return fig


def fig_bandwidth_vs_efficiency(theme, rows, rows_std):
    fig, ax = plt.subplots(figsize=(6.0, 4.4))
    ax.axvspan(MIN_EFFICIENCY * 100, 100, ymin=0, ymax=1, color=theme["series"][2], alpha=0.06)
    for i, sub in enumerate(CATALOGUE):
        rs = [r for r in rows if r["substrate"] == sub.name]
        ax.plot([r["efficiency"] * 100 for r in rs], [r["bandwidth_pct"] for r in rs],
                color=theme["series"][i], lw=1.6, label=sub.name)
        pts = [r for r in rows_std if r["substrate"] == sub.name]
        ax.scatter([r["efficiency"] * 100 for r in pts], [r["bandwidth_pct"] for r in pts], s=22,
                   color=theme["series"][i], edgecolor=theme["surface"], linewidth=0.8, zorder=3)
    ax.axhline(BW_REQUIRED * 100, color=theme["secondary"], ls=":", lw=1.2)
    ax.axvline(MIN_EFFICIENCY * 100, color=theme["secondary"], ls=":", lw=1.2)
    ax.text(81, BW_REQUIRED * 100 + 0.1, "spec: both lines", color=theme["secondary"], fontsize=8)
    ax.set_xlabel("Radiation efficiency (%)")
    ax.set_ylabel("VSWR<=2 bandwidth (%)")
    ax.set_title("No catalogue board reaches the top-right corner")
    ax.set_xlim(0, 100)
    ax.legend(frameon=False, fontsize=7.5, loc="upper left")
    fig.tight_layout()
    return fig


def fig_q_breakdown(theme):
    hs = np.linspace(0.2e-3, H_MAX_LAMBDA0 * LAMBDA0, 60)
    ds = [synthesize(F0, RT_DUROID_5880, h) for h in hs]
    fig, ax = plt.subplots(figsize=(5.8, 3.9))
    for i, (name, vals) in enumerate([
            ("$Q_{rad}$", [d.q_rad for d in ds]), ("$Q_c$ (copper)", [d.q_c for d in ds]),
            ("$Q_d$ (tan d)", [d.q_d for d in ds]), ("$Q_t$ (total)", [d.q_total for d in ds])]):
        ax.plot(hs * 1e3, vals, color=theme["series"][i], lw=2.2 if i == 3 else 1.5,
                ls="-" if i != 2 else "--", label=name)
    ax.set_yscale("log")
    ax.set_xlabel("Substrate thickness h (mm)")
    ax.set_ylabel("Quality factor")
    ax.set_title("Where the Q comes from (RT/duroid 5880, 2.44 GHz)")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    return fig


def fig_s11(theme, candidates):
    f = np.linspace(2.25e9, 2.65e9, 801)
    fig, ax = plt.subplots(figsize=(6.0, 3.9))
    ax.axvspan(BAND_LO / 1e9, BAND_HI / 1e9, color=theme["series"][2], alpha=0.09)
    ax.text((BAND_LO + BAND_HI) / 2e9, -1.2, "ISM band", ha="center", color=theme["secondary"], fontsize=8)
    for i, (d, label) in enumerate(candidates):
        s = 20 * np.log10(np.maximum(np.abs(d.s11(f)), 1e-6))
        ax.plot(f / 1e9, s, color=theme["series"][i], lw=1.8, label=label)
    vswr2 = 20 * math.log10(1 / 3)
    ax.axhline(vswr2, color=theme["secondary"], ls=":", lw=1.2)
    ax.text(2.26, vswr2 - 1.6, "VSWR = 2 (-9.54 dB)", color=theme["secondary"], fontsize=8)
    ax.set_ylim(-35, 0)
    ax.set_xlabel("Frequency (GHz)")
    ax.set_ylabel("$|S_{11}|$ (dB)")
    ax.set_title("Best catalogue board on each laminate")
    ax.legend(frameon=False, fontsize=7.5, loc="lower left")
    fig.tight_layout()
    return fig


def fig_patterns(theme, d):
    fig, axes = plt.subplots(1, 2, subplot_kw={"projection": "polar"}, figsize=(8.4, 4.0))
    k0 = d.k0
    ang = np.linspace(-np.pi / 2, np.pi / 2, 721)
    e = e_plane_pattern(ang, k0 * d.h, k0 * d.l_eff)
    h = h_plane_pattern(ang + np.pi / 2, k0 * d.h, k0 * d.w)
    for ax, pat, title in [(axes[0], e, "E-plane"), (axes[1], h, "H-plane")]:
        db = 20 * np.log10(np.maximum(pat / pat.max(), 10 ** (-30 / 20)))
        ax.plot(ang, db, color=theme["series"][0], lw=1.8)
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        ax.set_thetamin(-90)
        ax.set_thetamax(90)
        ax.set_rlim(-30, 0)
        ax.set_rticks([-30, -20, -10, 0])
        ax.set_title(title, fontsize=10, pad=8)
        ax.tick_params(labelsize=7.5, colors=theme["secondary"])
    fig.suptitle(f"Cavity-model patterns, {d.substrate.name}, h = {d.h*1e3:.3f} mm (dB, normalised)",
                 fontsize=9.5)
    fig.tight_layout()
    return fig


def write_csv(path, rows):
    keys = list(rows[0].keys())
    for r in rows[1:]:
        keys += [k for k in r if k not in keys]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main():
    os.makedirs("results", exist_ok=True)  # absent on a fresh CI checkout
    print(f"Target: {F0/1e9:.5f} GHz, VSWR<=2 bandwidth >= {BW_REQUIRED*100:.3f} %, "
          f"efficiency >= {MIN_EFFICIENCY:.0%}")
    rows = sweep(F0, H_SWEEP)
    rows_std = standard_board_options(F0, BW_REQUIRED, MIN_EFFICIENCY)
    ideal = solve_spec(F0, BW_REQUIRED, MIN_EFFICIENCY)
    write_csv("results/thickness_sweep.csv", rows)
    write_csv("results/catalogue_boards.csv", rows_std)
    write_csv("results/ideal_thickness_per_laminate.csv", ideal)

    best = {}
    for r in rows_std:  # best catalogue board per laminate: highest bandwidth among efficient ones,
        cur = best.get(r["substrate"])   # else highest bandwidth overall
        key = (r["meets_efficiency"], r["bandwidth_pct"])
        if cur is None or key > (cur["meets_efficiency"], cur["bandwidth_pct"]):
            best[r["substrate"]] = r
    subs = {s.name: s for s in CATALOGUE}
    candidates = [(synthesize(F0, subs[n], r["h_mm"] * 1e-3),
                   f"{n}, {r['h_mm']:.3f} mm ({r['efficiency']:.0%} eff.)")
                  for n, r in best.items() if n in (RT_DUROID_5880.name, RO4003C.name, FR4.name)]

    save_light_dark(lambda th: fig_bandwidth_efficiency_vs_h(th, rows, rows_std), "bandwidth_efficiency_vs_h")
    save_light_dark(lambda th: fig_bandwidth_vs_efficiency(th, rows, rows_std), "bandwidth_vs_efficiency")
    save_light_dark(fig_q_breakdown, "q_breakdown")
    save_light_dark(lambda th: fig_s11(th, candidates), "s11_best_boards")
    save_light_dark(lambda th: fig_patterns(th, synthesize(F0, RT_DUROID_5880, 3.175e-3)), "radiation_patterns")
    print("Wrote figures and CSVs to results/\n")

    print(f"{'laminate':20s} {'h (mm)':>7s} {'BW %':>6s} {'eff %':>6s} {'gain dBi':>8s} {'area mm2':>9s}  verdict")
    for r in rows_std:
        verdict = ("MEETS SPEC" if r["meets_spec"] else
                   "bandwidth only (loss)" if r["meets_bandwidth"] else "")
        if r["meets_bandwidth"] or r["h_mm"] == max(subs[r["substrate"]].standard_thicknesses_mm):
            print(f"{r['substrate']:20s} {r['h_mm']:7.3f} {r['bandwidth_pct']:6.2f} "
                  f"{r['efficiency']*100:6.1f} {r['gain_dBi']:8.2f} {r['area_mm2']:9.0f}  {verdict}")


if __name__ == "__main__":
    main()
