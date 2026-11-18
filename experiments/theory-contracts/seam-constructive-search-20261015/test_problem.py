"""Independent signed-charge, rational-clock and CPU/CUDA search-filter tests."""
from fractions import Fraction
import itertools
import unittest
from unittest.mock import patch

import numpy as np
import problem

try:
    import torch
except ImportError:
    torch = None


def product_oracle(n, c, m, d, chi, zeta):
    """Insert ordered integer generators, then reduce the final X square."""
    phase = c * sum(x*y for x, y in zip(chi, m))
    prefix = 0
    for nj, mj in zip(n, m):
        phase += 2 * nj * prefix
        prefix += mj
    charge = [x+y for x, y in zip(n, m)]
    if c and d:
        phase += zeta
        for j, value in enumerate(charge):
            phase += 2*j*value
        charge = [value+1 for value in charge]
    return charge, (c+d) % 2, phase


def associator_oracle(n, c, m, d, l, f, chi, zeta):
    nm, cd, p = product_oracle(n, c, m, d, chi, zeta)
    left, lg, q = product_oracle(nm, cd, l, f, chi, zeta)
    ml, df, r = product_oracle(m, d, l, f, chi, zeta)
    right, rg, s = product_oracle(n, c, ml, df, chi, zeta)
    return left, lg, right, rg, (p+q-r-s) % 4


def time_oracle(n, grade, alpha, erase_quarter=False, cyclic_charge=False):
    p = [Fraction(x) + Fraction(grade, 2) for x in n]
    image = p[4:] + p[:4]
    image = [x+int(j<4) for j, x in enumerate(image)]
    delta = [Fraction(j, 8) + alpha/8 - Fraction(1, 2) for j in range(8)]
    energy = lambda v: sum(x*x/2+d*x for x, d in zip(v, delta))
    q = sum(p)
    target_charge = q % 4 if cyclic_charge else q
    ward = energy(image)-energy(p)-(target_charge/2+Fraction(3,4)+alpha/2)
    # The fixed fine clock multiplies the j-weighted charge difference by1/8.
    character = sum(Fraction(j, 8)*(x-y) for j,(x,y) in enumerate(zip(image,p)))
    character += 0 if erase_quarter else alpha/2
    target = Fraction(6,8)+alpha/2+Fraction(target_charge % 2,2)
    return ward, sum(image)-q-4, (character-target) % 1


def candidates():
    positives = problem.exhaustive_generator_classification()
    mutants = []
    for g in positives:
        for j in range(8):
            m = {**g, "chi_exp": list(g["chi_exp"])}
            m["chi_exp"][j] = (m["chi_exp"][j]+1) % 4
            mutants.append(m)
        mutants.append({**g, "kappa_exp": (g["kappa_exp"]+1) % 4})
    mutants.append({**positives[0], "chi_exp": [0]*8})
    return positives, mutants


class SourceOracleTests(unittest.TestCase):
    def test_signed_normal_form_agrees_with_generator_insertion(self):
        rng = np.random.default_rng(8123)
        for _ in range(500):
            n, m = rng.integers(-20,21,(2,8))
            chi = rng.integers(0,4,8)
            zeta, c, d = (int(rng.integers(4)),int(rng.integers(2)),int(rng.integers(2)))
            charge, grade, phase = problem._normal_product(n,c,m,d,chi,zeta)
            expected = product_oracle(n.tolist(),c,m.tolist(),d,chi.tolist(),zeta)
            self.assertEqual((charge.tolist(),grade,phase),expected)

    def test_all_positive_classes_signed_associativity(self):
        rng = np.random.default_rng(442)
        for g in problem.exhaustive_generator_classification():
            for bits in itertools.product((0,1),repeat=3):
                for _ in range(20):
                    n,m,l = rng.integers(-100,101,(3,8)).tolist()
                    left,lg,right,rg,residual = associator_oracle(n,bits[0],m,bits[1],l,bits[2],g["chi_exp"],g["zeta_exp"])
                    self.assertEqual((left,lg,residual),(right,rg,0))

    def test_rational_time_and_character_signed_charges(self):
        rng = np.random.default_rng(917)
        self.assertEqual(problem._time_fixtures().shape,(36,3))
        self.assertTrue(np.array_equal(problem._time_fixtures(),np.zeros((36,3),dtype=np.int64)))
        for alpha in (Fraction(1,4),Fraction(3,4)):
            for grade in (0,1):
                for _ in range(100):
                    n = rng.integers(-100,101,8).tolist()
                    self.assertEqual(time_oracle(n,grade,alpha),(0,0,0))

    def test_source_invalid_time_mutants_are_detected_independently(self):
        n = [3,-2,7,1,-1,4,2,9]
        for alpha in (Fraction(1,4),Fraction(3,4)):
            self.assertNotEqual(time_oracle(n,1,alpha,erase_quarter=True)[2],0)
            self.assertNotEqual(time_oracle(n,1,alpha,cyclic_charge=True)[0],0)

    def test_random_sampler_is_not_a_calibration_seed(self):
        initial = problem.sample(np.random.default_rng(731),128)
        self.assertTrue(all(problem.valid(g) for g in initial))
        self.assertFalse(any(len(set(g["chi_exp"]))==1 and g["chi_exp"][0]%2 and g["kappa_exp"]%2==0 for g in initial))
        self.assertGreater(len({tuple(g["chi_exp"]) for g in initial}),120)

    def test_exhaustive_generator_constraints_have_exactly_four_classes(self):
        survivors = []
        for chi in itertools.product(range(4),repeat=8):
            if any((2*x-2)%4 for x in chi) or sum(chi)%4:
                continue
            if any((chi[i]-chi[j])%4 or (chi[i]+chi[j]-2)%4
                   for i in range(8) for j in range(i+1,8)):
                continue
            if any((chi[j]-chi[(j+4)%8])%4 for j in range(8)) or (-sum(chi[4:]))%4:
                continue
            for kappa in range(4):
                if 2*kappa%4==0:
                    survivors.append((chi,kappa))
        expected = [(tuple(g["chi_exp"]),g["kappa_exp"])
                    for g in problem.exhaustive_generator_classification()]
        self.assertEqual(survivors,expected)


@unittest.skipIf(torch is None,"Torch runtime unavailable")
class TensorEvaluatorTests(unittest.TestCase):
    def test_four_classes_pass_and_source_phase_mutants_fail(self):
        positives,mutants = candidates()
        rows = problem.evaluate(positives+mutants,"cpu")
        self.assertEqual(len(positives),4)
        self.assertTrue(all(r["feasible"] and r["objectives"]==[0]*4 for r in rows[:4]))
        self.assertTrue(all(not r["feasible"] for r in rows[4:]))
        self.assertTrue(all(r["details"]["associator_probes"]==200 and r["details"]["time_probes"]==36 for r in rows))

    def test_nonzero_time_fixture_is_consumed(self):
        old = dict(problem._CACHE)
        try:
            problem._CACHE.clear()
            with patch.object(problem,"_time_fixtures",return_value=np.array([[1,0,0]],dtype=np.int64)):
                rows = problem.evaluate(problem.exhaustive_generator_classification(),"cpu")
            self.assertTrue(all(not r["feasible"] and r["objectives"][3]==1 for r in rows))
        finally:
            problem._CACHE.clear()
            problem._CACHE.update(old)

    @unittest.skipUnless(torch is not None and torch.cuda.is_available(),"CUDA unavailable")
    def test_cpu_cuda_exact_parity(self):
        positives,mutants = candidates()
        gs = positives+mutants+problem.sample(np.random.default_rng(7144),128)
        cpu,cuda = problem.evaluate(gs,"cpu"),problem.evaluate(gs,"cuda:0")
        torch.cuda.synchronize()
        self.assertEqual([(r["objectives"],r["feasible"],r["descriptor"]) for r in cpu],
                         [(r["objectives"],r["feasible"],r["descriptor"]) for r in cuda])
        self.assertTrue(all(r["details"]["tensor_device"]=="cuda:0" for r in cuda))


if __name__ == "__main__":
    unittest.main()
