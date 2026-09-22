# Methodology: every equation, and where it comes from

Equation numbers in parentheses refer to C. A. Balanis, *Antenna Theory:
Analysis and Design*, 3rd ed., Wiley, 2005, Chapter 14. Everything here is
the transmission-line / cavity model of a rectangular patch in its dominant
TM010 mode. No full-wave electromagnetic solver (HFSS, CST, openEMS, FDTD)
is used anywhere in this repository.

## 1. Synthesis (`src/patch.py`, Balanis 14.2.1)

Width for efficient radiation (14-6):

```
W = c / (2 f0) * sqrt(2 / (eps_r + 1))
```

Effective permittivity of the microstrip formed by the patch, valid for
W/h > 1 (14-1):

```
eps_eff = (eps_r + 1)/2 + (eps_r - 1)/2 * (1 + 12 h/W)^(-1/2)
```

Hammerstad's fringing-field extension of each radiating edge (14-2):

```
dL = 0.412 h * (eps_eff + 0.3)(W/h + 0.264) / [(eps_eff - 0.258)(W/h + 0.8)]
```

Physical length that resonates at f0 (14-7):

```
L = c / (2 f0 sqrt(eps_eff)) - 2 dL,        L_eff = L + 2 dL
```

## 2. Radiating slots and input resistance (Balanis 14.2.2)

Each radiating edge is modelled as a slot of width W. With k0 = 2 pi / lambda0:

```
I1  = int_0^pi [ sin((k0 W/2) cos t) / cos t ]^2 sin^3 t dt           (14-53)
G1  = I1 / (120 pi^2)                                                  (14-12)
G12 = 1/(120 pi^2) int_0^pi [ sin((k0 W/2) cos t) / cos t ]^2
                             J0(k0 L sin t) sin^3 t dt                 (14-18a)
R_in(edge) = 1 / [ 2 (G1 + G12) ]      (dominant, odd-symmetry mode)   (14-17)
R_in(y0)   = R_in(edge) cos^2(pi y0 / L)                               (14-20a)
```

The integrand has a removable 0/0 at t = pi/2; the code writes
`sin(a cos t)/cos t` as `a * sinc(a cos t / pi)` so the limit is evaluated
exactly rather than by avoiding the point.

An inset feed can only move the feed point toward the centre of the patch,
where the resistance is lower. If the edge resistance is already below
Z0 = 50 ohm (which happens on very thin, lossy boards once loss is
included, see Section 4), no inset depth can match it: the design is fed at
the edge and flagged `matchable_by_inset = False`.

## 3. Directivity (Balanis 14.2.4)

Single slot (14-52), and two slots via the array factor (14-56):

```
D0   = (k0 W)^2 / I1
D_AF = 2 / (1 + g12),   g12 = G12 / G1
D2   = D0 * D_AF
```

Two slots from the full pattern integral (14-55), which is the expression
used for gain in this repository:

```
I2 = int_0^pi int_0^pi [ sin((k0 W/2) cos t) / cos t ]^2 sin^3 t
                        cos^2((k0 L_eff/2) sin t sin p) dt dp
D2 = (k0 W)^2 pi / I2
```

`verify.py` reproduces Balanis' Example 14.3 values for I1, D0, g12, D_AF
and D2 (14-56) to four significant figures. For I2, adaptive quadrature and
a 2000 x 2000 grid agree with each other on 3.56550 to twelve digits; the
book prints 3.59801, which corresponds to L_e = 1.057 cm rather than the
stated 1.068 cm. The resulting D2 differs from the book by 0.9 %, and the
check is set at 1 % with that explanation attached.

## 4. Loss, Q and bandwidth (derived here)

Balanis gives R_in only for the lossless patch. Bandwidth and efficiency
need the quality factors, derived below from the same TM010 cavity fields
so that they are consistent with the R_in above.

**Radiation Q.** Inside the cavity E_z = E0 cos(pi x / L). With peak
phasors, the time-averaged electric energy is

```
W_e = (eps / 4) int |E|^2 dV = (eps0 eps_r / 4) E0^2 h (W L / 2)
```

and at resonance W_m = W_e, so the stored energy is W = 2 W_e =
eps0 eps_r E0^2 h W L / 4. Each radiating edge carries voltage V0 = E0 h,
and the power radiated by the two coupled slots is the same power that
defines R_in in (14-17):

```
P_rad = V0^2 (G1 + G12)          [so that P_rad = V0^2 / (2 R_in)]
Q_rad = omega W / P_rad = omega eps0 eps_r W L / [ 4 h (G1 + G12) ]
```

**Conductor and dielectric Q** (14-88a, 14-88b):

```
Q_c = h sqrt(pi f mu0 sigma)          (sigma = 5.8e7 S/m, copper)
Q_d = 1 / tan(delta)
```

**Total Q, efficiency and edge resistance with loss:**

```
1/Q_t = 1/Q_rad + 1/Q_c + 1/Q_d        (surface-wave loss not included)
e_r   = Q_t / Q_rad                    (radiation efficiency)
R_edge(lossy) = R_in(edge) * e_r       (loss conductances in parallel)
```

**Frequency response.** Near resonance the patch at its feed behaves as a
parallel RLC circuit:

```
Z(f) = R / [ 1 + j Q_t (f/f0 - f0/f) ],     S11 = (Z - Z0) / (Z + Z0)
```

with R = Z0 for an inset-matched design and R = R_edge otherwise. The
VSWR <= 2 bandwidth is read numerically from where |S11| crosses 1/3. For a
matched design this reproduces the closed form

```
BW = (S - 1) / (Q_t sqrt(S)),   S = 2
```

and `verify.py` checks that the two agree.

**Why loss "helps" bandwidth.** Q_t falls whenever any loss term grows,
so a lossy laminate widens the VSWR bandwidth. That extra bandwidth is
dissipated in the board rather than radiated, which is why the design
study constrains efficiency as well as bandwidth.

## 5. The design question (`src/design_space.py`)

Target: the 2.4 GHz ISM band, 2400 to 2483.5 MHz, designed at its centre
(2.44175 GHz). Required VSWR <= 2 fractional bandwidth: 83.5 / 2441.75 =
3.42 %. Required radiation efficiency: >= 80 %.

For each laminate the study reports:

1. the thinnest board, as a continuous variable, that meets the bandwidth
   (a bracketing root-find, `brentq`, on bandwidth(h) - required; bandwidth
   rises monotonically with h once h is thick enough for radiation rather
   than loss to set Q_t, which holds throughout the region where the
   crossing occurs), within the 0.003 to 0.05 lambda0 range that Balanis
   gives as the usual thickness window for patches; and
2. every catalogue thickness a supplier actually stocks, flagged against
   both specifications.

The second is the realistic question, because a designer cannot order a
3.23 mm laminate.
