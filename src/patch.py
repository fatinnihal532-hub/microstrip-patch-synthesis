"""
Rectangular microstrip patch antenna: closed-form synthesis and analysis.

Everything in this module is the transmission-line / cavity model of a
rectangular patch operating in its dominant TM010 mode, as presented in

  C. A. Balanis, "Antenna Theory: Analysis and Design," 3rd ed., Wiley,
  2005, Chapter 14 (equation numbers below refer to that edition).

It is an analytical model, not a full-wave electromagnetic solve (no
HFSS, CST or FDTD is run anywhere in this repository). verify.py checks
this implementation number-for-number against Balanis' own worked
Examples 14.1, 14.2 and 14.3.

Quantities that the book gives only for the lossless case (radiation
Q, conductor Q, dielectric Q, and the resulting efficiency and
bandwidth) are derived here from the same cavity-model fields -- the
derivation is written out in docs/methodology.md.
"""
import math
from dataclasses import dataclass

import numpy as np
from scipy import integrate, optimize, special

from .substrates import Substrate

C0 = 299_792_458.0            # speed of light, m/s
MU0 = 4e-7 * math.pi          # vacuum permeability, H/m
EPS0 = 1.0 / (MU0 * C0 ** 2)  # vacuum permittivity, F/m
ETA_FACTOR = 120 * math.pi ** 2  # the 120*pi^2 that appears in (14-12), (14-18a)
SIGMA_COPPER = 5.8e7          # S/m, annealed copper


# ---------------------------------------------------------------------------
# 1. Synthesis: W and L for a target resonant frequency (Balanis 14.2.1)
# ---------------------------------------------------------------------------

def design_width(f0, eps_r):
    """(14-6): the width that gives efficient radiation."""
    return C0 / (2 * f0) * math.sqrt(2 / (eps_r + 1))


def effective_permittivity(eps_r, h, w):
    """(14-1), valid for W/h > 1."""
    if w / h <= 1:
        raise ValueError(f"W/h = {w / h:.3f} <= 1: outside the validity of (14-1)")
    return (eps_r + 1) / 2 + (eps_r - 1) / 2 * (1 + 12 * h / w) ** -0.5


def length_extension(eps_eff, h, w):
    """(14-2), Hammerstad's fringing-field length extension, per edge."""
    return 0.412 * h * ((eps_eff + 0.3) * (w / h + 0.264)) / ((eps_eff - 0.258) * (w / h + 0.8))


def design_length(f0, eps_eff, delta_l):
    """(14-7): physical length that resonates at f0 once both fringing
    extensions are accounted for."""
    return C0 / (2 * f0 * math.sqrt(eps_eff)) - 2 * delta_l


# ---------------------------------------------------------------------------
# 2. Radiating-slot conductances and input resistance (Balanis 14.2.2)
# ---------------------------------------------------------------------------

def _slot_pattern_sq(theta, k0w):
    """[sin((k0 W/2) cos t) / cos t]^2, written through sinc so that the
    removable 0/0 at theta = pi/2 is evaluated exactly (limit (k0W/2)^2)."""
    a = k0w / 2
    return (a * np.sinc(a * np.cos(theta) / np.pi)) ** 2


def slot_integral_i1(k0w):
    """(14-53) I1 = int_0^pi [sin((k0W/2)cos t)/cos t]^2 sin^3 t dt.
    The same integral appears in G1 (14-12) and in D0 (14-52)."""
    val, _ = integrate.quad(lambda t: _slot_pattern_sq(t, k0w) * math.sin(t) ** 3,
                            0, math.pi, limit=200)
    return val


def self_conductance(k0w):
    """(14-12): G1 = I1 / (120 pi^2)."""
    return slot_integral_i1(k0w) / ETA_FACTOR


def mutual_conductance(k0w, k0l):
    """(14-18a): G12 between the two radiating slots, a distance L apart."""
    val, _ = integrate.quad(
        lambda t: _slot_pattern_sq(t, k0w) * special.j0(k0l * math.sin(t)) * math.sin(t) ** 3,
        0, math.pi, limit=200)
    return val / ETA_FACTOR


def inset_depth_for(r_edge, r_target, length):
    """Invert (14-20a) Rin(y0) = Rin(0) cos^2(pi y0 / L) for y0."""
    if r_target > r_edge:
        raise ValueError(f"cannot match {r_target} ohm with an inset: edge resistance is only {r_edge:.1f} ohm")
    return length / math.pi * math.acos(math.sqrt(r_target / r_edge))


# ---------------------------------------------------------------------------
# 3. Directivity (Balanis 14.2.4)
# ---------------------------------------------------------------------------

def directivity_single_slot(k0w):
    """(14-52): D0 = (k0 W)^2 / I1."""
    return k0w ** 2 / slot_integral_i1(k0w)


def directivity_two_slot_af(k0w, g12_norm):
    """(14-56): D2 = D0 * D_AF, with D_AF = 2 / (1 + g12), g12 = G12/G1."""
    return directivity_single_slot(k0w) * 2 / (1 + g12_norm)


def slot_integral_i2(k0w, k0le):
    """(14-55a): I2 = int int [sin((k0W/2)cos t)/cos t]^2 sin^3 t
    cos^2((k0 Le/2) sin t sin p) dt dp over [0, pi] x [0, pi]."""
    val, _ = integrate.dblquad(
        lambda t, p: _slot_pattern_sq(t, k0w) * math.sin(t) ** 3
        * math.cos(k0le / 2 * math.sin(t) * math.sin(p)) ** 2,
        0, math.pi, 0, math.pi)
    return val


def directivity_two_slot(k0w, k0le):
    """(14-55): D2 = (k0 W)^2 pi / I2 -- the two-slot directivity from the
    full pattern integral, rather than the array-factor approximation."""
    return k0w ** 2 * math.pi / slot_integral_i2(k0w, k0le)


# ---------------------------------------------------------------------------
# 4. Cavity-model patterns (Balanis 14-43, 14-44), normalised to peak
# ---------------------------------------------------------------------------

def e_plane_pattern(phi, k0h, k0le):
    """(14-43): E-plane (theta = 90 deg) field magnitude vs phi, where phi
    is measured from broadside toward the radiating edges."""
    x = k0h / 2 * np.cos(phi)
    return np.abs(np.sinc(x / np.pi) * np.cos(k0le / 2 * np.sin(phi)))


def h_plane_pattern(theta, k0h, k0w):
    """(14-44): H-plane (phi = 90 deg) field magnitude vs theta."""
    x = k0h / 2 * np.sin(theta)
    z = k0w / 2 * np.cos(theta)
    return np.abs(np.sin(theta) * np.sinc(x / np.pi) * np.sinc(z / np.pi))


# ---------------------------------------------------------------------------
# 5. Loss, quality factor, bandwidth (derived; see docs/methodology.md)
# ---------------------------------------------------------------------------

def q_radiation(f0, eps_r, w, l, h, g_rad_total):
    """Q_rad = omega * W_stored / P_rad for the TM010 cavity:
    W_stored = eps E0^2 h W L / 4 (peak phasors, We = Wm at resonance),
    P_rad    = (E0 h)^2 (G1 + G12)  -- the same power that gives
               Rin = 1 / (2 (G1 + G12)) in (14-17).
    =>  Q_rad = omega eps0 eps_r W L / (4 h (G1 + G12)).
    """
    omega = 2 * math.pi * f0
    return omega * EPS0 * eps_r * w * l / (4 * h * g_rad_total)


def q_conductor(f0, h, sigma=SIGMA_COPPER):
    """(14-88a): Q_c = h sqrt(pi f mu0 sigma) for the two cavity walls."""
    return h * math.sqrt(math.pi * f0 * MU0 * sigma)


def q_dielectric(tan_delta):
    """(14-88b): Q_d = 1 / tan(delta)."""
    return 1.0 / tan_delta


def vswr_bandwidth(q_total, vswr=2.0):
    """Fractional bandwidth over which a parallel-RLC resonator matched at
    its centre stays below the given VSWR: (S - 1) / (Q_t sqrt(S))."""
    return (vswr - 1) / (q_total * math.sqrt(vswr))


# ---------------------------------------------------------------------------
# 6. The complete design object
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class PatchDesign:
    f0: float
    substrate: Substrate
    h: float
    w: float
    l: float
    eps_eff: float
    delta_l: float
    g1: float
    g12: float
    q_rad: float
    q_c: float
    q_d: float
    z0: float

    # -- derived quantities ----------------------------------------------
    @property
    def lam0(self):
        return C0 / self.f0

    @property
    def k0(self):
        return 2 * math.pi / self.lam0

    @property
    def l_eff(self):
        return self.l + 2 * self.delta_l

    @property
    def area_m2(self):
        return self.w * self.l

    @property
    def r_edge_lossless(self):
        """(14-17) with the + sign for the dominant (odd, TM010) mode."""
        return 1.0 / (2 * (self.g1 + self.g12))

    @property
    def q_total(self):
        return 1.0 / (1 / self.q_rad + 1 / self.q_c + 1 / self.q_d)

    @property
    def radiation_efficiency(self):
        """Fraction of accepted power that is radiated (conductor and
        dielectric loss only; surface-wave excitation is not modelled)."""
        return self.q_total / self.q_rad

    @property
    def r_edge(self):
        """Edge resistance including loss: the loss conductances appear in
        parallel with the radiation conductance, so R shrinks by the same
        factor Q_t/Q_rad as the efficiency."""
        return self.r_edge_lossless * self.radiation_efficiency

    @property
    def matchable_by_inset(self):
        """An inset can only move the feed toward the centre, where the
        resistance is lower, so it can match z0 only if R_edge >= z0. On a
        thin, lossy board R_edge (with loss) can fall below 50 ohm."""
        return self.r_edge >= self.z0

    @property
    def inset_depth(self):
        """Inset depth for a z0 match, or NaN if the patch has to be fed
        at the edge because R_edge < z0."""
        return inset_depth_for(self.r_edge, self.z0, self.l) if self.matchable_by_inset else float("nan")

    @property
    def feed_resistance(self):
        """Resonant resistance at the feed actually used: z0 through an
        inset if possible, otherwise the (too low) edge resistance."""
        return self.z0 if self.matchable_by_inset else self.r_edge

    @property
    def bandwidth_vswr2(self):
        """VSWR <= 2 fractional bandwidth, read numerically from S11 of
        the feed actually used (0 if even the centre frequency is worse
        than VSWR 2). For an inset-matched design this reproduces the
        closed form (S-1)/(Q_t sqrt(S)) -- verify.py checks that."""
        return self.numerical_bandwidth(2.0)

    def numerical_bandwidth(self, vswr):
        gamma_max = (vswr - 1) / (vswr + 1)
        excess = lambda f: abs(self.s11(f)) - gamma_max
        if excess(self.f0) > 0:
            return 0.0
        # |S11| grows monotonically away from f0 for a single parallel-RLC
        # resonance; bracket each band edge out to +/- 50 % of f0.
        f_lo = optimize.brentq(excess, 0.5 * self.f0, self.f0, xtol=1e-6 * self.f0)
        f_hi = optimize.brentq(excess, self.f0, 1.5 * self.f0, xtol=1e-6 * self.f0)
        return (f_hi - f_lo) / self.f0

    @property
    def directivity_af(self):
        """(14-56): array-factor approximation, D0 * 2/(1+g12)."""
        return directivity_two_slot_af(self.k0 * self.w, self.g12 / self.g1)

    @property
    def directivity(self):
        """(14-55): two-slot directivity from the full pattern integral,
        the more complete of Balanis' two expressions; used for gain."""
        return directivity_two_slot(self.k0 * self.w, self.k0 * self.l_eff)

    @property
    def gain(self):
        return self.radiation_efficiency * self.directivity

    # -- frequency response ------------------------------------------------
    def input_impedance(self, f):
        """Parallel-RLC model of the resonance seen at the feed:
        Z(f) = R / (1 + j Q_t (f/f0 - f0/f)), with R = feed_resistance.
        Feed-probe/line inductance is not included."""
        f = np.asarray(f, dtype=float)
        return self.feed_resistance / (1 + 1j * self.q_total * (f / self.f0 - self.f0 / f))

    def s11(self, f):
        z = self.input_impedance(f)
        return (z - self.z0) / (z + self.z0)


def synthesize(f0, substrate: Substrate, h, z0=50.0, w=None, sigma=SIGMA_COPPER):
    """Design a rectangular patch for resonance at f0 on the given
    substrate, following Balanis' design procedure (Sec. 14.2.1), then
    evaluate its conductances and Q factors. If w is given, it overrides
    the (14-6) "efficient radiator" width."""
    if h <= 0:
        raise ValueError("substrate thickness must be positive")
    w = design_width(f0, substrate.eps_r) if w is None else w
    eps_eff = effective_permittivity(substrate.eps_r, h, w)
    dl = length_extension(eps_eff, h, w)
    l = design_length(f0, eps_eff, dl)
    k0 = 2 * math.pi * f0 / C0
    g1 = self_conductance(k0 * w)
    g12 = mutual_conductance(k0 * w, k0 * l)
    qr = q_radiation(f0, substrate.eps_r, w, l, h, g1 + g12)
    return PatchDesign(f0=f0, substrate=substrate, h=h, w=w, l=l, eps_eff=eps_eff,
                       delta_l=dl, g1=g1, g12=g12, q_rad=qr,
                       q_c=q_conductor(f0, h, sigma), q_d=q_dielectric(substrate.tan_delta),
                       z0=z0)
