import unittest

from src.design_space import (sweep, solve_spec, standard_board_options, pareto_front,
                              minimum_thickness_for_bandwidth)
from src.patch import synthesize
from src.substrates import RT_DUROID_5880, CATALOGUE

F0 = 2.44175e9
BW = 83.5e6 / F0


class TestDesignSpace(unittest.TestCase):
    def test_sweep_row_count(self):
        rows = sweep(F0, [1.0e-3, 2.0e-3, 3.0e-3])
        self.assertEqual(len(rows), 3 * len(CATALOGUE))

    def test_minimum_thickness_meets_spec_exactly(self):
        h = minimum_thickness_for_bandwidth(F0, RT_DUROID_5880, BW, 0.4e-3, 6e-3)
        self.assertAlmostEqual(synthesize(F0, RT_DUROID_5880, h).bandwidth_vswr2 / BW, 1.0, places=4)

    def test_unreachable_bandwidth_returns_none(self):
        self.assertIsNone(minimum_thickness_for_bandwidth(F0, RT_DUROID_5880, 0.5, 0.4e-3, 6e-3))

    def test_solve_spec_flags_lossy_fr4(self):
        rows = {r["substrate"]: r for r in solve_spec(F0, BW, 0.80)}
        self.assertTrue(rows["RT/duroid 5880"]["meets_spec"])
        self.assertFalse(rows["FR-4"]["meets_spec"])

    def test_catalogue_boards_have_consistent_flags(self):
        for r in standard_board_options(F0, BW, 0.80):
            self.assertEqual(r["meets_spec"], r["meets_bandwidth"] and r["meets_efficiency"])

    def test_pareto_front_hand_checked_example(self):
        rows = [{"g": 7, "a": 10}, {"g": 6, "a": 5}, {"g": 5, "a": 8}, {"g": 4, "a": 2}]
        # (5, 8) is beaten by (6, 5): more gain and less area. The other three
        # each win on one objective.
        front = pareto_front(rows, "g", "a")
        self.assertEqual(sorted((r["g"], r["a"]) for r in front), [(4, 2), (6, 5), (7, 10)])


if __name__ == "__main__":
    unittest.main()
