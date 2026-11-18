#!/usr/bin/env python3
"""Exact conditional Hardy/CAR diagnostic; no derivation of P1 is claimed."""
from __future__ import annotations
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations_with_replacement, product
import json
from math import comb, factorial
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
MANDATORY_REPO_FILES = (
    'tfpt_1_architecture_e8.tex',
    '_archive/tfpt-45/source_extracts/01_boundary_kernel_source.tex',
    'experiments/lean4-carrier-rigidity/TfptCarrier/SeamWindingInterface.lean',
)
MANDATORY_LOCAL_FILES = ('search.py','PROOF.txt','lean/SourceReconstruction.lean','lean-receipt.json')
ROOTS = tuple(map(F, ('-1/2', '-1/4', '0', '1/4', '1/2')))
z, x = sp.symbols('z x', real=True)


def canon(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=False)


def digest(value):
    return sha256(canon(value).encode()).hexdigest()


def rational(q):
    q = F(q)
    return sp.Rational(q.numerator, q.denominator)


def reduced(numerator, denominator=()):
    a, b = Counter(map(F, numerator)), Counter(map(F, denominator))
    if any(abs(t) >= 1 for t in a.keys() | b.keys()):
        raise ValueError('Blaschke roots must be real and strictly inside the disk')
    common = a & b
    a.subtract(common)
    b.subtract(common)
    return tuple(sorted(a.elements())), tuple(sorted(b.elements()))


def bfactor(a):
    a = rational(a)
    return (z-a)/(1-a*z)


def symbol(a, b):
    return sp.cancel(sp.prod(bfactor(t) for t in a)/sp.prod(bfactor(t) for t in b))


def negative_parts(a, b):
    """Principal parts of inside poles, which are precisely negative Laurent modes."""
    u = symbol(a, b)
    parts = []
    for root, multiplicity in sorted(Counter(b).items()):
        t = rational(root)
        regular = sp.cancel(u*(z-t)**multiplicity)
        for k in range(1, multiplicity+1):
            degree = multiplicity-k
            coefficient = F(sp.diff(regular, z, degree).subs(z, t)/factorial(degree))
            if coefficient:
                parts.append((root, k, coefficient))
    return parts


def kernel(a, k, b, ell):
    """Exact L2 inner product of (z-a)^-k and (z-b)^-ell on |z|=1."""
    p, q = k-1, ell-1
    return sum((F(factorial(p+q-j), factorial(j)*factorial(p-j)*factorial(q-j))
                * b**(p-j)*a**(q-j)/(1-a*b)**(p+q-j+1)
                for j in range(min(p,q)+1)), F(0))


def leakage(parts):
    return sum((c*d*kernel(a,k,b,ell)
                for a,k,c in parts for b,ell,d in parts), F(0))


def negative_coefficient(parts, n):
    if n < 1:
        raise ValueError('n indexes coefficient of z^-n and must be positive')
    return sum((c*comb(n-1,k-1)*a**(n-k)
                for a,k,c in parts if n >= k), F(0))


def phase_derivative(a, b):
    def poisson(t):
        t = rational(t)
        return (1-t*t)/(1-2*t*x+t*t)
    return sp.cancel(sum(map(poisson,a), sp.S.Zero)-sum(map(poisson,b), sp.S.Zero))


def derivative_data(a, b):
    def bounds(t):
        lo = (1-abs(t))/(1+abs(t))
        return lo, 1/lo
    lower = sum((bounds(t)[0] for t in a),F(0))-sum((bounds(t)[1] for t in b),F(0))
    upper = sum((bounds(t)[1] for t in a),F(0))-sum((bounds(t)[0] for t in b),F(0))
    derivative = phase_derivative(a,b)
    at_zero = F(derivative.subs(x,1))
    at_pi = F(derivative.subs(x,-1))
    return {'formula_in_cos_theta':str(derivative),
            'certified_lower_bound':str(lower), 'certified_upper_bound':str(upper),
            'bounds_sharp_by_endpoints': min(at_zero,at_pi)==lower and max(at_zero,at_pi)==upper,
            'at_theta_zero':str(at_zero), 'at_theta_pi':str(at_pi),
            'strict_orientation_preserving_certified':lower>0,
            'orientation_reversal_witness':at_zero<0 or at_pi<0,
            'original_D_increment_identically_one':sp.cancel(derivative-1)==0}


def same_source_energy(a,b):
    """Quasifree state of U with the original D_alpha fixed, never conjugated."""
    signed=[(1,t) for t in a]+[(-1,t) for t in b]
    winding=len(a)-len(b)
    mean_square=sum((s*t*(1+u*v)/(1-u*v)
                     for s,u in signed for t,v in signed),F(0))
    variance=mean_square-winding*winding
    if variance<0:
        raise ArithmeticError('Negative exact phase-derivative variance')
    return {'original_B_fixed':True,
            'state':'Bogoliubov transform of the same negative sea by scalar multiplication U',
            'premises':'given Hardy/CAR source, sea state and original B=N+alpha; smooth implementable U',
            'relative_charge_expectation':winding,
            'phase_derivative_mean_square':str(mean_square),
            'phase_derivative_variance':str(variance),
            'energy_expectation_by_alpha':{
                str(alpha):str(mean_square/2+(alpha-F(1,2))*winding)
                for alpha in (F(1,4),F(3,4))},
            'energy_excess_above_winding_one_floor':str(variance/2) if winding==1 else None,
            'attains_winding_one_energy_minimum':winding==1 and variance==0}


def evaluate(numerator, denominator=()):
    a,b = reduced(numerator,denominator)
    parts = negative_parts(a,b)
    leak = leakage(parts)
    # This is an algebraic pole test, not a finite Fourier cutoff test.
    inner = not b
    if inner != (not parts) or (leak == 0) != inner:
        raise ArithmeticError('Reduced pole classification and exact norm disagree')
    winding = len(a)-len(b)
    data = {'numerator':list(map(str,a)), 'denominator':list(map(str,b)),
            'winding':winding, 'negative_principal_parts':[
                {'pole':str(t),'order':k,'coefficient':str(c)} for t,k,c in parts],
            'negative_coefficients_first_six': [str(negative_coefficient(parts,n)) for n in range(1,7)],
            'vacuum_negative_leakage_norm_squared':str(leak),
            'compression_identity':'(C U C)* (C U C) = C - C U* (1-C) U C on ran(C)',
            'hardy_isometry_in_assumed_model':inner,
            'pure_defect_one_in_assumed_model':inner and winding==1,
            'same_source_energy':same_source_energy(a,b),
            'original_B':derivative_data(a,b)}
    data['id'] = digest({'numerator':data['numerator'],'denominator':data['denominator']})
    return data


def enumerate_symbols():
    """Bounded diagnostics only: winding 1, at most 3 uncancelled factors."""
    seen = set()
    for p,q in ((1,0),(2,1)):
        for a in combinations_with_replacement(ROOTS,p):
            for b in combinations_with_replacement(ROOTS,q):
                pair = reduced(a,b)
                if pair not in seen:
                    seen.add(pair)
                    yield pair


def reconstructed_operator(a):
    """Check B_u on its rational eigenvectors; normalization cancels here."""
    t = rational(a)
    b = bfactor(a)
    k = 1/(1-t*z)
    logarithmic_k = t*z/(1-t*z)
    poisson = sp.cancel(z*sp.diff(b,z)/b)
    checks = []
    for n in range(-2,3):
        g = k*b**n
        checks.append(sp.cancel((z*sp.diff(g,z)-logarithmic_k*g)/poisson-n*g)==0)
    return {'root':str(a), 'F_action':'F_a f = sqrt(1-a^2)/(1-a*z) * f(b_a(z))',
            'F_preserves_Hardy_projection':'analytic weighted composition, conditional on supplied Hardy space',
            'B_u_action':'B_u = P_a^-1 (z*d/dz - a*z/(1-a*z)) + alpha',
            'eigenvector_modes_checked':list(range(-2,3)),
            'eigenvector_identities':checks,
            'defect_vector':'k_a(z)=sqrt(1-a^2)/(1-a*z)',
            'defect_vector_norm_squared':str((1-a*a)/(1-a*a)),
            'Hardy_reproducing_kernel_zero':str(sp.cancel(b.subs(z,t))),
            'defect_orthogonality_reason':'<b_a z^n,k_a> = sqrt(1-a^2) b_a(a) a^n = 0, n>=0',
            'B_u_equals_supplied_D_alpha':a==0,
            'alpha_values':['1/4','3/4'],
            'quarter_monodromy_exponents':['1/4','3/4'],
            'full_unitarity_and_completeness':'analytic proof obligation; finite modes do not establish it'}


def occupancy(modes, bits):
    return {n for n,bit in zip(modes,bits) if n>=0 and bit}, {n for n,bit in zip(modes,bits) if n<0 and bit}


def invariants(particles, holes, alpha):
    if any(n<0 for n in particles) or any(n>=0 for n in holes):
        raise ValueError('Particles must be nonnegative and holes strictly negative')
    return len(particles)-len(holes), sum((F(n)+alpha for n in particles),F(0))-sum((F(n)+alpha for n in holes),F(0))


def shift_sea(particles, holes, steps=1):
    """Translate the occupied infinite sea, retaining finite deviations from n<0."""
    lo = min([-1,*particles,*holes])+min(steps,0)-2
    hi = max([0,*particles,*holes])+max(steps,0)+2
    occupied = lambda n: (n<0 and n not in holes) or n in particles
    new_particles = {n for n in range(0,max(1,hi+1)) if occupied(n-steps)}
    new_holes = {n for n in range(min(lo,-1),0) if not occupied(n-steps)}
    return new_particles,new_holes


def field_action(particles, holes, mode, create=False):
    """Exact sea-basis CAR action, ordered from larger modes to smaller modes."""
    occupied = (mode<0 and mode not in holes) or mode in particles
    if occupied == create:
        return None
    above = max(0,-mode-1)+sum(n>mode for n in particles)-sum(n>mode for n in holes)
    p,h = set(particles),set(holes)
    if mode>=0:
        (p.add if create else p.remove)(mode)
    else:
        (h.remove if create else h.add)(mode)
    return (-1)**above,p,h


def field_audit():
    records=[]
    for alpha in (F(1,4),F(3,4)):
        for bits in product((0,1),repeat=6):
            p,h=occupancy((-3,-2,-1,0,1,2),bits)
            q,energy=invariants(p,h,alpha)
            for mode in (-3,-2,-1,0,1,2):
                for create in (False,True):
                    action=field_action(p,h,mode,create)
                    shifted=field_action(*shift_sea(p,h),mode+1,create)
                    if action is None:
                        if shifted is not None: raise ArithmeticError('Field seam intertwiner failed')
                        continue
                    sign,p1,h1=action
                    if shifted!=(sign,*shift_sea(p1,h1)):
                        raise ArithmeticError('Fermionic seam sign/intertwiner failed')
                    q1,e1=invariants(p1,h1,alpha)
                    delta=1 if create else -1
                    if q1-q!=delta or e1-energy!=delta*(mode+alpha):
                        raise ArithmeticError('Field charge/time commutator failed')
                    records.append([str(alpha),list(bits),mode,create,sign,q1,str(e1)])
    return {'nonzero_field_arrows_checked':len(records),'arrows_digest':digest(records),
            'CAR_representation':'assumed sea CAR, exact occupation signs',
            'intertwiner':'R c_n R* = c_(n+1)',
            'time_commutator':'[H_alpha,c_n]=-(n+alpha)c_n',
            'charge_commutator':'[Q,c_n]=-c_n'}


def core_audit():
    rows=[]
    modes=(-3,-2,-1,0,1,2)
    for alpha in (F(1,4),F(3,4)):
        for bits in product((0,1),repeat=len(modes)):
            particles,holes=occupancy(modes,bits)
            q,h=invariants(particles,holes,alpha)
            p1,h1=shift_sea(particles,holes)
            q1,e1=invariants(p1,h1,alpha)
            if q1!=q+1 or e1!=h+q+alpha or (h-alpha*q).denominator!=1:
                raise ArithmeticError('Infinite-sea field/state/time identity failed')
            # Inverse shift catches spurious finite-window/cyclic implementations.
            if shift_sea(p1,h1,-1)!=(particles,holes):
                raise ArithmeticError('Sea translation failed inverse test')
            rows.append({'alpha':str(alpha),'particles':sorted(particles),'holes':sorted(holes),
                         'Q':q,'H':str(h),'shifted_Q':q1,'shifted_H':str(e1)})
    vacuum=[]
    for alpha in (F(1,4),F(3,4)):
        for k in range(-4,5):
            p,h=shift_sea(set(),set(),k)
            charge,energy=invariants(p,h,alpha)
            expected=F(k*(k-1),2)+alpha*k
            if charge!=k or energy!=expected:
                raise ArithmeticError('Charged vacuum orbit failed')
            vacuum.append({'alpha':str(alpha),'k':k,'charge':charge,'energy':str(energy)})
    return {'premise':'pure rank-one Hardy shift, minimal unitary dilation, actual B=N+alpha, charged CAR representation',
            'field_map':'beta_U(c_n)=c_(n+1)', 'vacuum_occupation':'n<0',
            'configuration_count':len(rows),'configurations_per_alpha':64,
            'rows_digest':digest(rows), 'rows':rows,
            'field_action_audit':field_audit(),
            'charge_identity':'R* Q R = Q+1', 'time_identity':'R* H_alpha R = H_alpha+Q+alpha',
            'quarter_identity':'exp(2*pi*i*H_alpha)=exp(2*pi*i*alpha*Q)',
            'vacuum_orbit':vacuum, 'seam_preserves_vacuum':False,
            'U_order':'infinite','alpha_transition_order':4,
            'full_finite_core_statement':'conditional analytic proof reference; bounded audit is not the proof'}


def source_obligations():
    # Evidence is absent, not inferred from operator checks in an assumed model.
    premises = {
        'Hardy_projection_C_from_P1':None,
        'actual_multiplication_U_from_P1':None,
        'actual_C_U_leakage_zero_from_P1':None,
        'actual_shift_pure_and_minimal_dilation_from_P1':None,
        'original_B_equals_reconstructed_B':None,
        'CAR_and_field_identification_from_P1':None,
        'P2_carrier_to_actual_field_dictionary':None,
        'physical_vacuum_state_identification_from_P1':None,
        'physical_time_generator_identification_from_P1':None}
    premises['minimum_energy_state_selection_from_P1']=None
    unmet=[key for key,evidence in premises.items() if not evidence]
    return {'evidence':premises,'unmet':unmet,
            'route_specific_obligations':{
                'Wold_isometry_route':['actual_C_U_leakage_zero_from_P1','actual_shift_pure_and_minimal_dilation_from_P1','original_B_equals_reconstructed_B'],
                'same_source_variational_route':['physical_time_generator_identification_from_P1','CAR_and_field_identification_from_P1',
                    'physical_vacuum_state_identification_from_P1','minimum_energy_state_selection_from_P1']},
            'outcome':'NO_P1_DERIVATION' if unmet else 'REQUIRES_EXTERNAL_PROOF_REVIEW',
            'reason':'Hardy/CAR and source operators are assumed inputs of this diagnostic, not consequences of P1.'}


def build():
    rows=[evaluate(a,b) for a,b in enumerate_symbols()]
    controls={'canonical_z':evaluate([0]),'single_b_half':evaluate([F(1,2)]),
              'oriented_noninner_quarter':evaluate([0,0],[F(1,4)]),
              'folding_noninner_half':evaluate([0,0],[F(1,2)])}
    return {'schema':'conditional-source-reconstruction-v1',
        'domain':{'roots':list(map(str,ROOTS)),'max_total_factors':3,'winding':1,
                  'global_phase':'fixed to 1; phase does not affect classification',
                  'nonzero_leakage_scope':'Excludes Wold-isometry route only. Every smooth unimodular symbol here still defines a conditional Bogoliubov quasifree CAR state.',
                  'bounded_enumeration_is_full_family_proof':False},
        'classification':{'diagnostic_count':len(rows),
             'hardy_isometries':sum(r['hardy_isometry_in_assumed_model'] for r in rows),
             'preserve_original_D_increment':sum(r['original_B']['original_D_increment_identically_one'] for r in rows),
             'same_source_energy_minima':sum(r['same_source_energy']['attains_winding_one_energy_minimum'] for r in rows),
             'full_family_conditional_statement':'For continuous scalar unimodular symbols, winding-one Hardy isometries are single Blaschke factors. Keeping original D_alpha covariance selects a=0, up to constant phase.',
             'proof_reference':'PROOF.txt (analytic argument; not certified by this finite enumeration)'},
        'same_source_variational_selection':{
            'given_source':'original C, B=N+alpha, charged CAR and negative sea remain fixed',
            'energy_formula':'<H_alpha>=1/2 mean(p^2)+(alpha-1/2)w',
            'winding_one_formula':'<H_alpha>=alpha+1/2 Var(p)',
            'conditional_full_family_minimizer':'For smooth scalar H1 winding-one symbols, energy equality iff p=1, hence u=constant phase times z.',
            'selection_from_P1':'unproved; minimum-energy selection is an additional principle',
            'proof_reference':'PROOF.txt; finite factor enumeration is only a diagnostic'},
        'candidates':rows,'controls':controls,
        'reconstructed_operator_controls':[reconstructed_operator(t) for t in (F(0),F(1,4),F(1,2))],
        'conditional_field_state_time':core_audit(), 'source_obligations':source_obligations()}


def package(repo_root=None, local_root=None):
    repo_root=Path(repo_root) if repo_root is not None else REPO
    local_root=Path(local_root) if local_root is not None else HERE
    data=build()
    # Mandatory originals bind the actual source audit. No fallback to old hashes.
    paths=[repo_root/name for name in MANDATORY_REPO_FILES]
    paths += [local_root/name for name in MANDATORY_LOCAL_FILES]
    snapshot=local_root/'atlas_snapshot.json'
    if snapshot.exists(): paths.append(snapshot)
    sources={path.relative_to(repo_root).as_posix():sha256(path.read_bytes()).hexdigest()
             for path in paths}
    receipt=json.loads((local_root/'lean-receipt.json').read_text())
    lean_hash=sha256((local_root/'lean/SourceReconstruction.lean').read_bytes()).hexdigest()
    if receipt['source']['sha256']!=lean_hash:
        raise ValueError('Lean receipt source hash disagrees with the actual Lean source')
    return {'sources':sources,'payload':data,'payload_sha256':digest(data)}


def verify(document):
    expected=package()
    if document!=expected:
        raise ValueError('Replay mismatch: result mutation, changed source, or changed analytic proof reference')
    return True


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=('search','verify'))
    parser.add_argument('--output',type=Path,default=HERE/'results.json')
    args=parser.parse_args()
    if args.command=='search':
        document=package()
        args.output.write_text(json.dumps(document,indent=2,sort_keys=True)+'\n')
        print(json.dumps({'outcome':document['payload']['source_obligations']['outcome'],
                          'counts':document['payload']['classification'],
                          'payload_sha256':document['payload_sha256']}))
    else:
        verify(json.loads(args.output.read_text()))
        print('VERIFIED exact replay; conditional diagnostic, NO_P1_DERIVATION')

if __name__=='__main__':
    main()
