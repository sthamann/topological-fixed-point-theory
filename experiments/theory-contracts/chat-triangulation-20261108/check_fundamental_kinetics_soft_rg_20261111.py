"""Three source-typed tests: kinetic symmetry, soft universality and gauge matching."""

import json
from pathlib import Path
import sys

import numpy as np
import sympy as sp


root = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(root / 'verification'))
import v1002_spin10_mirror_projector as carrier
import v246_unification_data_crosscheck as gauge_data
import v249_ps_unification as ps
import v252_full_finite_triple as matter
import v310_carrier_sm_anomaly as sm


def multiplet(occupied):
    color = len(occupied & {0, 1, 2})
    weak = len(occupied & {3, 4})
    return {(0, 0): 'nuc', (2, 0): 'uc', (1, 1): 'Q', (0, 2): 'ec',
            (2, 2): 'dc', (3, 1): 'L'}[color, weak]


def edge_matrix(names, edges):
    basis = {name: sp.eye(len(names)).row(i) for i, name in enumerate(names)}
    return sp.Matrix.vstack(*(basis[a] - basis[b] for a, b in sorted(edges)))


def kinetics():
    names = [name for name, _, _ in sm.GEN]
    half_spinor = [i for i in range(32) if i.bit_count() % 2 == 0]
    labels = [multiplet({j for j in range(5) if state & (16 >> j)}) for state in half_spinor]
    assert {name: labels.count(name) for name in names} == {name: size for name, size, _ in sm.GEN}
    gamma = carrier.gamma_matrices_10d()
    d5 = set()
    connecting_generator = None
    quarks = {'Q', 'uc', 'dc'}
    leptons = set(names) - quarks
    for a in range(10):
        for b in range(a+1, 10):
            generator = (gamma[a] @ gamma[b])[np.ix_(half_spinor, half_spinor)]
            for i, j in zip(*np.nonzero(generator), strict=True):
                first, second = labels[i], labels[j]
                if first == second:
                    continue
                assert abs(generator[i, j]) == 1
                d5.add(tuple(sorted((first, second))))
                if first in quarks and second in leptons and connecting_generator is None:
                    connecting_generator = {'clifford_indices_zero_based': [a, b], 'between': [second, first]}

    # Use the actual native operator support, with all allowed Yukawa blocks active.
    unit, zero = np.eye(3), np.zeros((3, 3))
    def fields(kind, chirality):
        return [slot for (name, side, _, _), slot in matter.IDX.items()
                if (name, side) == (kind, chirality)]

    yukawa = matter.full_D(unit, unit, unit, unit, zero)
    edges = set()
    for kind, left, right in [('u', 'Q', 'uc'), ('d', 'Q', 'dc'),
                              ('e', 'L', 'ec'), ('nu', 'L', 'nuc')]:
        assert np.max(abs(yukawa[np.ix_(fields(kind, 'L'), fields(kind, 'R'))])) == 1
        edges.add(tuple(sorted((left, right))))
    sm_constraints = edge_matrix(names, edges)
    assert len(sm_constraints.nullspace()) == 2
    assert edge_matrix(names, d5).rank() == 5
    assert connecting_generator is not None
    assert edge_matrix(names, edges | {tuple(sorted(connecting_generator['between']))}).rank() == 5

    # Quark/lepton speeds are an allowed kinetic deformation of the supplied SM representation.
    labels_48 = np.zeros(48)
    for (name, _, _, _), slot in matter.IDX.items():
        labels_48[slot] = 2 if name in ('u', 'd') else 1
    diagonal = np.tile(labels_48, 2)
    v = np.diag(diagonal)
    native = [matter.full_D(unit, 2*unit, 3*unit, 4*unit, unit), matter.Kswap(),
              matter.grading(), matter.rep(1j, np.array([[0, 1], [-1, 0]]),
                                           np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]]))]
    assert all(np.array_equal(v @ operator, operator @ v) for operator in native)

    # Even sign flips of the five weight axes are Weyl(D5) operations.
    # In the spin weights they act by even bit flips, with the spin lift's
    # phases cancelling on this diagonal kinetic coefficient.
    vq, vl = sp.symbols('v_q v_l')
    values = dict(zip(half_spinor, (vq if label in quarks else vl for label in labels), strict=True))
    averaged = {sp.expand(sum(values[state ^ flip] for flip in half_spinor)/16) for state in half_spinor}
    assert averaged == {3*vq/4 + vl/4}
    return {'SM_yukawa_speed_dimension': 2, 'D5_compatible_speed_dimension_per_family': 1,
            'one_sufficient_native_D5_generator': connecting_generator,
            'D5_links': sorted(d5), 'native_96_quark_lepton_speed_split_commutators': 'zero',
            'even_D5_weight_flip_average_per_family': str(averaged.pop())}


def soft_universality():
    # Each line is an elastic transfer with a nonzero SM vertex connecting its legs.
    # On asymptotic states the leading longitudinal soft factor is (k_i-k_j)*transfer.
    sectors = ['electron', 'neutrino', 'quark', 'photon', 'gluon', 'Higgs', 'weak']
    edges = [('electron', 'neutrino'), ('electron', 'quark'), ('electron', 'photon'),
             ('quark', 'gluon'), ('electron', 'Higgs'), ('electron', 'weak')]
    constraints = edge_matrix(sectors, edges)
    assert constraints.nullspace() == [sp.ones(len(sectors), 1)]
    metric = sp.diag(1, -1, -1, -1)
    rest, light, exchange = (sp.Matrix(x) for x in ([1, 0, 0, 0], [1, 0, 0, 1], [0, 1, 0, 0]))
    massive_stress, light_stress = rest*rest.T, light*light.T
    assert massive_stress*metric*exchange == light_stress*metric*exchange == sp.zeros(4, 1)
    trace_mass = sp.trace(metric*massive_stress)
    trace_light = sp.trace(metric*light_stress)
    scalar_residue = trace_mass*trace_light
    tensor_residue = sp.trace(metric*massive_stress*metric*light_stress) - trace_mass*trace_light/2
    assert scalar_residue == 0 and tensor_residue == 1
    return {'sector_count': len(sectors), 'ward_rank': constraints.rank(),
            'remaining_overall_spin2_couplings': len(constraints.nullspace()), 'edges': edges,
            'normalized_mass_radiation_tensor_contraction': str(tensor_residue),
            'mass_radiation_trace_scalar_contraction': str(scalar_residue)}


def gauge_rank():
    b = sp.Matrix([sp.Rational(x.numerator, x.denominator) for x in sm.beta_coeffs(3)])
    a = sp.Matrix([sp.Rational(str(x)) for x in gauge_data.AINV])
    assert np.allclose(np.array(b, dtype=float).ravel(), ps.BSM)
    common = sp.ones(3, 1)
    sm_normal = common.cross(b)
    assert sm_normal.dot(a) != 0

    scalaron_a = a - b * sp.log(sp.Float(ps.M_S)/sp.Float(ps.M_Z))/(2*sp.pi)
    branches = {}
    for name, coefficients in [('minimal_16', ps.B_MIN), ('16_plus_45', ps.B_45)]:
        b4, b_l, b_r = (sp.Rational(float(x)).limit_denominator(100) for x in coefficients)
        ps_direction = sp.Matrix([(3*b_r + 2*b4)/5, b_l, b4])
        normal = common.cross(ps_direction)
        assert normal.dot(common) == 0 and normal.dot(ps_direction) == 0
        residual = float(normal.dot(scalaron_a) / normal.norm())
        assert abs(residual) > 1e-3
        natural_ps, natural_gut, _, ok = ps.ps_solve(coefficients)
        assert ok
        branches[name] = {'unfitted_threshold_over_scalaron': natural_ps / ps.M_S,
                          'original_GUT_GeV': natural_gut,
                          'unit_normal_inverse_coupling_residual_at_scalaron': residual,
                          'normal': list(map(str, normal))}

    f, g = sp.symbols('f g', positive=True)
    beta = [coefficient*g**3/(16*sp.pi**2) - f*g for coefficient in b]
    assert sp.sign(b[0]) == 1 and sp.sign(b[1]) == sp.sign(b[2]) == -1
    assert all(sp.diff(expr, g).subs(g, 0) == -f for expr in beta)
    assert sp.simplify(sp.diff(beta[0], g).subs(g, sp.sqrt(16*sp.pi**2*f/b[0])) - 2*f) == 0
    return {'native_beta': list(map(str, b)), 'SM_only_exact_normal_residual': str(sm_normal.dot(a)),
            'scalaron_GeV': ps.M_S, 'scalaron_threshold_branches': branches,
            'universal_antiscreening_positive_f': {'interacting_U1_slope': '+2 f',
                                                  'gaussian_SU2_SU3_slopes': '-f each'}}


print(json.dumps({'kinetics': kinetics(), 'soft_ward': soft_universality(), 'gauge_matching': gauge_rank()}, indent=2))
