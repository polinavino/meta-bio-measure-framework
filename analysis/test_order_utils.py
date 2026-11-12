"""Mathematical regression checks for the shared order/sampling implementation."""
import unittest
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from order_utils import average_rank, dominance, pair_disagreement, exact_rank_draws, _chain, _executable
from kinase_controls import candidate
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "inflammation/analysis"))
from sepsis_preprocessing import discovery_standardize


class OrderTests(unittest.TestCase):
    def test_validation_expression_does_not_fit_standardization(self):
        x = pd.DataFrame([[1., 2., 20.], [5., 7., -10.]], columns=["d1", "d2", "v1"])
        z, mean, sd = discovery_standardize(x, ["d1", "d2"])
        changed = x.copy(); changed["v1"] *= 1000
        zz, mm, ss = discovery_standardize(changed, ["d1", "d2"])
        pd.testing.assert_series_equal(mean, mm)
        pd.testing.assert_series_equal(sd, ss)
        pd.testing.assert_frame_equal(z[["d1", "d2"]], zz[["d1", "d2"]])

    def test_candidate_domain_and_scale_conditions(self):
        x = np.array([[5., 5.], [5.8, 5.], [7., 5.2]])
        result = candidate(x, 5., 6.)
        self.assertTrue(np.isnan(result[:2]).all())
        self.assertTrue(np.isfinite(result[2]))
        with self.assertRaises(ValueError):
            candidate(x, 5., 6., scale=.5)
        weak = np.column_stack([x[2:], [6.]])
        self.assertLessEqual(candidate(weak, 5., 6.)[0], result[2] + 1e-10)

    def test_ties_are_not_reversals(self):
        x = np.array([[0, 0], [1, 0], [0, 1], [0, 0]])
        i, j = np.triu_indices(len(x), 1)
        strict, ties = pair_disagreement(x, i, j)
        ge = dominance(x)
        np.testing.assert_array_equal(strict, ~(ge | ge.T)[i, j])
        self.assertEqual(strict.sum(), 1)
        self.assertTrue(ties.any())

    def test_exact_partial_order_and_equivalent_objects(self):
        x = np.array([[3, 2], [2, 3], [1, 1], [1, 1]])
        np.testing.assert_allclose(average_rank(x), [1.5, 1.5, 3, 3])

    def test_row_and_monotone_invariance(self):
        x = np.array([[1, 4], [2, 2], [4, 1], [5, 5]])
        a = average_rank(x)
        np.testing.assert_allclose(a, average_rank(x[::-1])[::-1])
        np.testing.assert_allclose(a, average_rank(np.column_stack([np.exp(x[:, 0]), x[:, 1] * 1000])))

    def test_native_kernel_against_enumerated_distribution(self):
        # Exactly three extensions, with the first object constrained before second.
        before = np.zeros((3, 3), bool)
        before[0, 1] = True
        exact = exact_rank_draws(before)
        draws = _chain(_executable(), before, [0, 1, 2], 100, 12000, 1, 77)
        self.assertTrue(np.all(draws[:, 0] < draws[:, 1]))
        np.testing.assert_allclose(draws.mean(0), exact.mean(0), atol=.035)
        states, counts = np.unique(draws, axis=0, return_counts=True)
        self.assertEqual(len(states), 3)
        np.testing.assert_allclose(counts / len(draws), np.full(3, 1/3), atol=.02)


if __name__ == "__main__":
    unittest.main()
