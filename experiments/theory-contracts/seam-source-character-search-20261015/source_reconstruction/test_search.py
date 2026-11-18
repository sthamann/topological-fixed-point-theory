"""Exact regression controls; no physical derivation is asserted."""
import copy
from fractions import Fraction as F
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
import sympy as sp

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('source_reconstruction_search',HERE/'search.py')
s=importlib.util.module_from_spec(spec)
spec.loader.exec_module(s)

class OperatorTests(unittest.TestCase):
    def test_cancellation_and_rejection(self):
        self.assertEqual(s.reduced([0,F(1,2)],[F(1,2)]),((F(0),),()))
        self.assertEqual(s.evaluate([0,F(1,2)],[F(1,2)]),s.evaluate([0]))
        for bad in ('1','-1','3/2'):
            with self.assertRaises(ValueError): s.evaluate([bad])

    def test_primary_oriented_counterexample(self):
        a=F(1,4)
        row=s.evaluate([0,0],[a])
        self.assertEqual(row['winding'],1)
        self.assertEqual(row['vacuum_negative_leakage_norm_squared'],'15/4096')
        self.assertFalse(row['hardy_isometry_in_assumed_model'])
        self.assertEqual(row['original_B']['certified_lower_bound'],'1/3')
        self.assertTrue(row['original_B']['strict_orientation_preserving_certified'])
        self.assertTrue(row['original_B']['bounds_sharp_by_endpoints'])
        # Independent rational long-division identity and geometric Laurent tail.
        z=s.z
        polynomial=-s.rational(a)*z*z+(1-s.rational(a)**2)*z+s.rational(a)*(1-s.rational(a)**2)
        residue=a*a*(1-a*a)
        self.assertEqual(sp.cancel(s.symbol((F(0),F(0)),(a,))-polynomial-s.rational(residue)/(z-s.rational(a))),0)
        self.assertEqual(residue*residue/(1-a*a),F(15,4096))
        self.assertEqual(s.negative_coefficient(s.negative_parts((F(0),F(0)),(a,)),4),residue*a**3)

    def test_folding_control(self):
        row=s.evaluate([0,0],[F(1,2)])
        self.assertEqual(row['vacuum_negative_leakage_norm_squared'],'3/64')
        self.assertEqual(row['original_B']['at_theta_zero'],'-1')
        self.assertTrue(row['original_B']['orientation_reversal_witness'])

    def test_repeated_inside_pole(self):
        # 1/b_0^2 = z^-2 has exactly unit leakage, no z^-1 component.
        parts=s.negative_parts((),(F(0),F(0)))
        self.assertEqual(s.leakage(parts),1)
        self.assertEqual(s.negative_coefficient(parts,1),0)
        self.assertEqual(s.negative_coefficient(parts,2),1)
        # 1/b_a is unimodular with constant -a, hence negative norm 1-a^2.
        a=F(1,4)
        self.assertEqual(s.leakage(s.negative_parts((),(a,))),1-a*a)
        self.assertEqual(s.leakage(s.negative_parts((),(a,a))),1-a**4)

    def test_bound_enumeration(self):
        rows=[s.evaluate(a,b) for a,b in s.enumerate_symbols()]
        self.assertEqual(len(rows),55)
        self.assertEqual(len({r['id'] for r in rows}),55)
        self.assertEqual(sum(r['hardy_isometry_in_assumed_model'] for r in rows),5)
        self.assertEqual(sum(r['original_B']['original_D_increment_identically_one'] for r in rows),1)
        self.assertEqual(sum(r['same_source_energy']['attains_winding_one_energy_minimum'] for r in rows),1)
        for row in rows:
            self.assertEqual(row['winding'],1)
            self.assertEqual(row['hardy_isometry_in_assumed_model'],F(row['vacuum_negative_leakage_norm_squared'])==0)

    def test_reconstructed_not_original_time(self):
        for a in (F(0),F(1,4),F(1,2)):
            row=s.reconstructed_operator(a)
            self.assertTrue(all(row['eigenvector_identities']))
            self.assertEqual(row['defect_vector_norm_squared'],'1')
            self.assertEqual(row['Hardy_reproducing_kernel_zero'],'0')
            self.assertEqual(row['B_u_equals_supplied_D_alpha'],a==0)
            self.assertEqual(row['quarter_monodromy_exponents'],['1/4','3/4'])
        self.assertTrue(s.evaluate([F(1,2)])['hardy_isometry_in_assumed_model'])
        self.assertFalse(s.evaluate([F(1,2)])['original_B']['original_D_increment_identically_one'])

    def test_infinite_sea_core(self):
        core=s.core_audit()
        self.assertEqual(core['configuration_count'],128)
        self.assertFalse(core['seam_preserves_vacuum'])
        for row in core['rows']:
            alpha=F(row['alpha'])
            self.assertEqual(row['shifted_Q'],row['Q']+1)
            self.assertEqual(F(row['shifted_H']),F(row['H'])+row['Q']+alpha)
            self.assertEqual((F(row['H'])-alpha*row['Q']).denominator,1)
        # Excitations far outside the audit window must translate as well.
        p,h={100},{-91}
        self.assertEqual(s.shift_sea(*s.shift_sea(p,h,7),-7),(p,h))
        self.assertEqual(s.shift_sea(set(),set()),({0},set()))
        self.assertEqual(s.shift_sea(set(),set(),-1),(set(),{-1}))

    def test_unmet_source_obligations(self):
        result=s.source_obligations()
        self.assertEqual(len(result['unmet']),10)
        self.assertIn('P2_carrier_to_actual_field_dictionary',result['unmet'])
        self.assertIn('minimum_energy_state_selection_from_P1',result['unmet'])
        self.assertEqual(result['outcome'],'NO_P1_DERIVATION')
        self.assertEqual(set(result['unmet']),set(k for k,v in result['evidence'].items() if v is None))

    def test_same_source_energy_from_independent_fourier_sums(self):
        # For coefficients c_l, sea translation gives energy weight
        # l(l-1)/2 + alpha*l. Geometric tails below are summed exactly.
        for a in (F(1,4),F(1,2),F(-1,4)):
            t=a*a
            sum_n=t/(1-t)**2
            sum_n2=t*(1+t)/(1-t)**3
            # b_a: c_0=-a, c_n=(1-a^2)a^(n-1), n>=1.
            q_inner=(1-t)**2*sum_n/t
            square_inner=(1-t)**2*sum_n2/t
            # z^2/b_a: c_2=-a, c_1=1-a^2, c_0=a(1-a^2),
            # c_-n=(1-a^2)a^(n+1), n>=1.
            q_noninner=2*t+(1-t)**2-(1-t)**2*t*sum_n
            square_noninner=4*t+(1-t)**2+(1-t)**2*t*sum_n2
            self.assertEqual(q_inner,1)
            self.assertEqual(q_noninner,1)
            for numerator,denominator,q,square in (
                ([a],[],q_inner,square_inner),([0,0],[a],q_noninner,square_noninner)):
                result=s.evaluate(numerator,denominator)['same_source_energy']
                self.assertTrue(result['original_B_fixed'])
                self.assertEqual(F(result['phase_derivative_mean_square']),square)
                self.assertEqual(F(result['energy_excess_above_winding_one_floor']),t/(1-t))
                for alpha in (F(1,4),F(3,4)):
                    expected=(square-q)/2+alpha*q
                    self.assertEqual(F(result['energy_expectation_by_alpha'][str(alpha)]),expected)
        primary=s.evaluate([0,0],[F(1,4)])['same_source_energy']
        self.assertEqual(primary['energy_expectation_by_alpha']['1/4'],'19/60')
        self.assertEqual(primary['energy_excess_above_winding_one_floor'],'1/15')

    def test_actual_fermionic_field_signs(self):
        # Direct CAR control on all six-mode occupations, with no cyclic boundary.
        from itertools import product
        def composed(p,h,m,n,first_create,second_create):
            first=s.field_action(p,h,n,first_create)
            if first is None: return None
            sign,p1,h1=first
            second=s.field_action(p1,h1,m,second_create)
            if second is None: return None
            other,p2,h2=second
            return sign*other,frozenset(p2),frozenset(h2)
        for bits in product((0,1),repeat=6):
            p,h=s.occupancy((-3,-2,-1,0,1,2),bits)
            for m,n in product(range(-3,3),repeat=2):
                for first_create in (False,True):
                    first=composed(p,h,m,n,first_create,False)
                    second=composed(p,h,n,m,False,first_create)
                    terms={}
                    for term in (first,second):
                        if term:
                            sign,p1,h1=term
                            key=(p1,h1)
                            terms[key]=terms.get(key,0)+sign
                    terms={k:v for k,v in terms.items() if v}
                    expected={(frozenset(p),frozenset(h)):1} if first_create and m==n else {}
                    self.assertEqual(terms,expected)
        self.assertEqual(s.field_audit()['nonzero_field_arrows_checked'],768)

    def test_fresh_process_replay_and_mutations(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/'results.json'
            command=[sys.executable,str(HERE/'search.py')]
            subprocess.run(command+['search','--output',str(out)],check=True,capture_output=True,text=True)
            subprocess.run(command+['verify','--output',str(out)],check=True,capture_output=True,text=True)
            original=json.loads(out.read_text())
            for field in ('payload','sources'):
                mutated=copy.deepcopy(original)
                if field=='payload':
                    mutated['payload']['controls']['oriented_noninner_quarter']['vacuum_negative_leakage_norm_squared']='0'
                else:
                    mutated['sources'][next(iter(mutated['sources']))]='0'*64
                out.write_text(json.dumps(mutated))
                process=subprocess.run(command+['verify','--output',str(out)],capture_output=True,text=True)
                self.assertNotEqual(process.returncode,0)
                self.assertIn('Replay mismatch',process.stderr)
            # Mutate the actual executable source, not merely its recorded hash.
            clone=Path(tmp)/'clone'
            local=clone/HERE.relative_to(s.REPO)
            for path in original['sources']:
                target=clone/path
                target.parent.mkdir(parents=True,exist_ok=True)
                shutil.copyfile(s.REPO/path,target)
            altered=local/'search.py'
            unchanged=altered.read_text()
            altered.write_text(unchanged+'\n# changed source\n')
            out.write_text(json.dumps(original))
            process=subprocess.run([sys.executable,str(altered),'verify','--output',str(out)],capture_output=True,text=True)
            self.assertNotEqual(process.returncode,0)
            self.assertIn('Replay mismatch',process.stderr)
            # The real original is never touched: mutate only its temporary clone.
            altered.write_text(unchanged)
            source=clone/s.MANDATORY_REPO_FILES[0]
            source.write_text(source.read_text()+'\n% changed original source\n')
            process=subprocess.run([sys.executable,str(altered),'verify','--output',str(out)],capture_output=True,text=True)
            self.assertNotEqual(process.returncode,0)
            self.assertIn('Replay mismatch',process.stderr)
            source.unlink()
            with self.assertRaises(FileNotFoundError):
                s.package(clone,local)

if __name__=='__main__': unittest.main()
