"""Exact topology and rank16 root controls; no invented genus5 theta values."""
import unittest
from f1_double_torus_genus5_parent import (
    overlay_genus, projection_ramification_count,
    castelnuovo_severi_trigonal_bound, rank16_roots, root_graph_components, controls,
)

A=(-4,-3,-2,-1)

class DoubleTorusParentTests(unittest.TestCase):
    def test_disjoint_branch_double_torus_is_genus_five(self):
        self.assertEqual(overlay_genus(A,(1,2,3,4)),5)
        self.assertEqual(projection_ramification_count(A,(1,2,3,4)),8)
        self.assertEqual(2*5-2,8)

    def test_shared_branch_switch_and_disconnected_equal_covers(self):
        self.assertEqual([overlay_genus(A, tuple(A[:k])+tuple(range(1,5-k)))
                          for k in range(4)], [5,4,3,2])
        self.assertIsNone(overlay_genus(A,A))

    def test_bielliptic_genus_five_excludes_degree_three_map(self):
        self.assertEqual(castelnuovo_severi_trigonal_bound(),4)
        self.assertGreater(overlay_genus(A,(1,2,3,4)),castelnuovo_severi_trigonal_bound())

    def test_two_rank16_parent_lattices_have_same_480_roots(self):
        e=rank16_roots('E8xE8')
        d=rank16_roots('D16plus')
        self.assertEqual(len(e),len(d))
        self.assertEqual(len(e),480)
        self.assertTrue(all(sum(x*x for x in r)==8 for r in e))
        self.assertTrue(all(sum(x*x for x in r)==8 for r in d))

    def test_same_root_number_but_different_pairing_graph(self):
        self.assertEqual(root_graph_components(rank16_roots('E8xE8')),(240,240))
        self.assertEqual(root_graph_components(rank16_roots('D16plus')),(480,))

    def test_parent_theta_separation_requires_literature_not_bare_enumeration(self):
        r=controls()
        self.assertFalse(r['published_genus5_partition_difference_evaluated_here'])
        self.assertFalse(r['original_Q0_found'])
        self.assertFalse(r['F1_solved'])

if __name__=='__main__':
    unittest.main()
