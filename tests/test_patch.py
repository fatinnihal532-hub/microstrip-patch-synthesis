import math
import unittest

import numpy as np

from src.patch import (synthesize, design_width, effective_permittivity, length_extension,
                       inset_depth_for, vswr_bandwidth, e_plane_pattern, h_plane_pattern, C0)
from src.substrates import Substrate, RT_DUROID_5880, FR4

F0 = 2.44175e9


class TestSynthesis(unittest.TestCase):
    def test_designed_patch_resonates_at_target(self):
        # (14-7) is built so that c / (2 Leff sqrt(eps_eff)) returns f0.
        d = synthesize(F0, RT_DUROID_5880, 1.575e-3)
        fr = C0 / (2 * d.l_eff * math.sqrt(d.eps_eff))
        self.assertAlmostEqual(fr / F0, 1.0, places=12)

    def test_eps_eff_lies_between_one_and_eps_r(self):
        for eps_r in (2.2, 4.4, 10.2):
            w = design_width(F0, eps_r)
            e = effective_permittivity(eps_r, 1.5e-3, w)
            self.assertTrue(1.0 < e < eps_r)

    def test_eps_eff_rejects_narrow_strip(self):
        with self.assertRaises(ValueError):
            effective_permittivity(2.2, h=2e-3, w=1e-3)

    def test_fringing_extension_is_order_of_h(self):
        # Hammerstad's dL is a fraction of h (0.3h-0.5h for common geometries).
        h = 1.575e-3
        w = design_width(F0, 2.2)
        dl = length_extension(effective_permittivity(2.2, h, w), h, w)
        self.assertTrue(0.3 * h < dl < 0.6 * h)

    def test_higher_permittivity_shrinks_the_patch(self):
        small = synthesize(F0, Substrate("hi-k", 10.2, 0.002), 1.27e-3)
        large = synthesize(F0, RT_DUROID_5880, 1.27e-3)
        self.assertLess(small.area_m2, large.area_m2)


class TestFeedAndBandwidth(unittest.TestCase):
    def test_inset_formula_round_trips(self):
        r_edge, l = 230.0, 0.04
        y0 = inset_depth_for(r_edge, 50.0, l)
        self.assertAlmostEqual(r_edge * math.cos(math.pi * y0 / l) ** 2, 50.0, places=9)

    def test_inset_cannot_raise_resistance(self):
        with self.assertRaises(ValueError):
            inset_depth_for(30.0, 50.0, 0.04)

    def test_thin_lossy_board_is_flagged_unmatchable(self):
        d = synthesize(F0, RT_DUROID_5880, 0.127e-3)
        self.assertFalse(d.matchable_by_inset)
        self.assertTrue(math.isnan(d.inset_depth))

    def test_bandwidth_formula_matches_numerical_s11(self):
        d = synthesize(F0, RT_DUROID_5880, 1.575e-3)
        self.assertAlmostEqual(d.bandwidth_vswr2 / vswr_bandwidth(d.q_total), 1.0, places=3)

    def test_efficiency_is_a_fraction(self):
        for sub in (RT_DUROID_5880, FR4):
            for h in (0.5e-3, 1.6e-3, 3.2e-3):
                e = synthesize(F0, sub, h).radiation_efficiency
                self.assertTrue(0.0 < e < 1.0)

    def test_lossier_laminate_is_less_efficient(self):
        low_loss = synthesize(F0, Substrate("A", 4.4, 0.001), 1.6e-3)
        high_loss = synthesize(F0, Substrate("B", 4.4, 0.02), 1.6e-3)
        self.assertGreater(low_loss.radiation_efficiency, high_loss.radiation_efficiency)


class TestPatterns(unittest.TestCase):
    def test_broadside_is_the_pattern_maximum(self):
        d = synthesize(F0, RT_DUROID_5880, 1.575e-3)
        ang = np.linspace(-np.pi / 2, np.pi / 2, 721)
        e = e_plane_pattern(ang, d.k0 * d.h, d.k0 * d.l_eff)
        h = h_plane_pattern(ang + np.pi / 2, d.k0 * d.h, d.k0 * d.w)
        self.assertEqual(int(np.argmax(e)), 360)
        self.assertEqual(int(np.argmax(h)), 360)

    def test_h_plane_nulls_at_grazing(self):
        # The sin(theta) factor in (14-44) forces a null along the ground plane.
        d = synthesize(F0, RT_DUROID_5880, 1.575e-3)
        self.assertLess(h_plane_pattern(np.array([0.0]), d.k0 * d.h, d.k0 * d.w)[0], 1e-12)


if __name__ == "__main__":
    unittest.main()
