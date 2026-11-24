"""Exact negative controls for potential clauses and their scope."""
from itertools import combinations
import json
import unittest
import sympy as sp

from process_metric_potential_certificates import (
    BETA_V, BETA_CUB, BETA_DET, SHIFT_MATRIX, metric_potential_audit,
    metric_potential_certificate, metric_potential_class_proof,
    metric_vertex_probe_matrix, star_quartic_orbits,
)


class MetricPotentialCertificateTests(unittest.TestCase):
    def test_six_candidates_are_classified_without_samples(self):
        audit = metric_potential_audit()
        rows = {x["name"]: set(x["clauses"]) for x in audit["candidates"]}
        self.assertEqual(len(rows["V"]), 6)
        self.assertEqual(rows["V + det(I-X)/2"], rows["V"]-{"quartic_as_V"})
        self.assertEqual(rows["V - det(I-X)/2"], rows["V"]-{"quartic_as_V"})
        cubic = {"flat_background", "fierz_pauli_mass", "massless_direction_free"}
        self.assertEqual(rows["V + cubic/2"], cubic)
        self.assertEqual(rows["V + cubic/2 + det(I-X)/2"], cubic)
        self.assertEqual(rows["beta1 only"], set())
        self.assertFalse(audit["samples_decide_clauses"])
        json.dumps(audit)

    def test_mass_combination_does_not_certify_a_nonflat_background(self):
        d = metric_potential_certificate((0, 1, 0, 0, 0))
        self.assertEqual(d["fierz_pauli_combination"], "1")
        self.assertFalse(d["clause_truths"]["fierz_pauli_mass"])

    def test_exact_input_rejects_rounded_or_malformed_betas(self):
        for beta in [(-2, 0.5, 0, 0.5, -2), (1, 2, 3), (True, 0, 0, 0, 0), (sp.oo, 0, 0, 0, 0)]:
            with self.assertRaises(ValueError):
                metric_potential_certificate(beta)

    def test_general_two_parameter_family_and_exchange(self):
        alpha, tau = sp.symbols("alpha tau", real=True)
        beta = [v+alpha*c+tau*d for v, c, d in zip(BETA_V, BETA_CUB, BETA_DET)]
        row = metric_potential_certificate(beta)
        self.assertTrue(row["flat_unit_mass_deformation"]["exact_decomposition"])
        self.assertEqual(row["flat_unit_mass_deformation"]["alpha_cubic"], "alpha")
        self.assertEqual(row["flat_unit_mass_deformation"]["tau_determinant"], "tau")
        undirected = metric_potential_certificate([v+tau*d for v, d in zip(BETA_V, BETA_DET)])
        self.assertTrue(undirected["clause_truths"]["undirected_edge"])
        self.assertTrue(undirected["clause_truths"]["cubic_as_V"])
        self.assertFalse(undirected["clause_truths"]["quartic_as_V"])

    def test_general_rank_proofs_and_higher_orders(self):
        self.assertTrue(metric_potential_class_proof()["all_checks_pass"])
        self.assertEqual(metric_vertex_probe_matrix(3).rank(), 4)
        for k in (4, 5, 9):
            self.assertEqual(metric_vertex_probe_matrix(k).rank(), 5)

    def test_quartic_tensor_monomials_from_independent_diagonal_expansion(self):
        eps = sp.Symbol("eps")
        ds = sp.symbols("d1:5")
        betas = sp.symbols("b0:5")
        roots = [sum(sp.binomial(sp.Rational(1, 2), k)*(eps*d)**k for k in range(5)) for d in ds]
        V = sum(betas[n]*sum((sp.prod(c) for c in combinations(roots, n)), sp.Integer(0)) for n in range(5))
        quartic = sp.Poly(sp.expand(V).coeff(eps, 4), *ds)
        monomials = [ds[0]**4, ds[0]**3*ds[1], ds[0]**2*ds[1]*ds[2], sp.prod(ds)]
        predicted = metric_vertex_probe_matrix(4)*sp.Matrix(betas)
        for r, monomial in enumerate(monomials, 1):
            self.assertEqual(sp.expand(quartic.coeff_monomial(monomial)-predicted[r]), 0)
        self.assertEqual(predicted[0], -sp.Rational(5, 128)*(SHIFT_MATRIX*sp.Matrix(betas))[0])

    def test_equal_metric_direction_can_vanish_without_flatness(self):
        d = metric_potential_certificate((1, 0, 0, 0, -1))
        self.assertTrue(d["clause_truths"]["massless_direction_free"])
        self.assertFalse(d["clause_truths"]["flat_background"])

    def test_tau_is_independent_per_depth_under_mu4(self):
        for m in (1, 2, 4):
            d = star_quartic_orbits(m)
            self.assertEqual(d["independent_tau_before_mu4"], 4*m)
            self.assertEqual(d["mu4_constraint_rank"], 3*m)
            self.assertEqual(d["independent_tau_after_mu4"], m)
            self.assertTrue(d["homogeneous_tau_across_depth_is_extra"])

    def test_canonical_quartic_normalization_and_action_sign(self):
        mu2, M, tau = sp.symbols("mu2 M tau", positive=True)
        ds = sp.symbols("d0:4")
        determinant = sp.prod(ds)
        density_quartic = tau*sp.prod(-2*d/M for d in ds)/16
        self.assertEqual(sp.simplify(-mu2*M*M*density_quartic+tau*mu2*determinant/(M*M)), 0)

    def test_external_exact_candidates_need_no_lane_import(self):
        d = metric_potential_audit({"custom": BETA_V})
        self.assertEqual(len(d["candidates"]), 1)
        self.assertTrue(d["candidates"][0]["all_orders_same_potential_as_V"])


if __name__ == "__main__":
    unittest.main()
