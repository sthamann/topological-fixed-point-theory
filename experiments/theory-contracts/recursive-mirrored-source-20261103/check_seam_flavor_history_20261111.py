"""Test historical Flavor formulas and candidate source readouts, without fitting."""
import json
from itertools import permutations, product

import numpy as np
import sympy as s

from tfpt_explorer.flavor_path_transport import _integrate, exact_clock_dictionary
from tfpt_constants import phi0
from v141_deck_selection import G, TA, P
from v94_sheet_diamond import Q, K, L, A, ONE, pi_left, pi_right


PAIRS = [(0, 1), (0, 2), (1, 2)]


def exterior(matrix):
    return s.Matrix(3, 3, lambda i, j: matrix[list(PAIRS[i]), list(PAIRS[j])].det())


def rep(perm):
    # H_1 = Z^4 / Z(1,1,1,1); fourth small loop is minus the other three.
    marks = [*s.eye(3).columnspace(), -s.ones(3, 1)]
    return s.Matrix.hstack(*(marks[perm[i]] for i in range(3)))


def monodromy(radius, rtol):
    mats = []
    for k in range(4):
        p = 1j**k
        endpoint = (1 - radius)*p
        leg = _integrate(lambda t: (t*endpoint, endpoint), rtol, rtol/100)[0]

        def arc(t):
            q = -p*radius*np.exp(2j*np.pi*t)
            return p+q, 2j*np.pi*q

        loop = _integrate(arc, rtol, rtol/100)[0]
        mats.append(np.linalg.solve(leg, loop @ leg))
    basis = [np.diag(row).astype(complex) for row in np.eye(3)]
    for i, j in PAIRS:
        for x in [1, 1j]:
            mat = np.zeros((3, 3), complex)
            mat[i, j], mat[j, i] = x, np.conj(x)
            basis.append(mat)
    U = np.diag([1, 1j, -1j])
    equations = np.stack([np.concatenate([(m.conj().T @ h @ m-h).ravel()
                                         for m in mats+[U]]) for h in basis], axis=1)
    _, sigma, vh = np.linalg.svd(np.vstack([equations.real, equations.imag]), full_matrices=False)
    H = sum(c*b for c, b in zip(vh[-1], basis))
    H /= np.trace(H).real/3
    I = np.eye(3)
    W = np.diag(np.sqrt(np.diag(H).real))
    actual = W @ mats[0] @ np.linalg.inv(W)
    target = np.array(exact_clock_dictionary()['M'], dtype=complex)
    assert np.linalg.eigvalsh(H)[0] > 0 and sigma[-2] > 0.1 and sigma[-1] < 1e-8
    assert np.linalg.norm(actual-target) < 1e-8
    return {'H_diagonal_trace_3': np.diag(H).real.tolist(), 'equation_svd_tail': sigma[-2:].tolist(),
            'unitary_defect': max(float(np.linalg.norm(m.conj().T @ H @ m-H)) for m in mats),
            'M0_difference': float(np.linalg.norm(actual-target)),
            'local_order3': max(float(np.linalg.norm(m@m@m-I)) for m in mats),
            'cyclic_product': float(np.linalg.norm(mats[0]@mats[1]@mats[2]@mats[3]-I))}


def test():
    rat = s.Rational
    targets = [rat(55, 117), rat(34, 47), rat(3, 26)]
    actual = [rat(5*sum(map(abs, pi_left(K))), 9*13),
              (sum(pi_left(L))-sum(pi_left(Q)))/s.Integer(1+sum(pi_left(L))),
              rat(3, sum(map(abs, pi_right(K))))]
    assert targets == actual
    assert [rat(16, 7), rat(4, 3), rat(7, 6)] == [rat(2)/(1-rat(1, 2)**3),
        rat(3, 2)/(1+rat(1, 2)**3), rat(32, 9)/(rat(16, 7)*rat(4, 3))]
    d = s.Symbol('d')
    cos = [1, rat(1, 2), -rat(1, 2), -1, -rat(1, 2), rat(1, 2)]
    D = lambda y, r: y*y+d*d-2*y*d*cos[r]
    gcds = {}
    for u, v in product([rat(1, 3), rat(2, 3), s.Integer(1)], repeat=2):
        f = [s.Poly(D(v, rd)-a*D(u, ru), d).monic()
             for a, ru, rd in zip(targets, (1, 3, 0), (1, 5, 2))]
        gcds[f'{u},{v}'] = str(s.gcd(s.gcd(f[0], f[1]), f[2]).as_expr())
    assert set(gcds.values()) == {'1'}
    ud_equation = D(rat(1, 3), 1)-targets[0]*D(rat(2, 3), 1)
    delta_roots = s.solve(ud_equation, d)
    assert len(delta_roots) == 2 and all(root.is_real for root in delta_roots)
    historical_table = [[.466, 1.085, .899], [.440, .646, .986], [.936, .943, .480]]
    old_products = [float(np.prod(sector)) for sector in historical_table]
    expected_volume_ratio = float(s.prod(targets))*(1-float(phi0))**2/float(phi0)
    assert all(p < .9 for p in old_products) and abs(expected_volume_ratio-1) > .3
    # If kinetic normalization is the only missing amplitude, the shared
    # Q_L volume cancels and the right Gram determinant ratio is fixed.
    target_right_gram = expected_volume_ratio**2
    delta_test = [3/5+float(phi0)/6, ((794-7*np.sqrt(9961))/2187)**(1/6), .5]
    mobius = [((1-z)/(1+z))**2 for z in delta_test]
    c6_lepton_e_mu = [4*(1+z*z-z)/(1+z*z+z) for z in delta_test]
    assert abs(c6_lepton_e_mu[2]-12/7) < 1e-15

    deck = rep((1, 2, 3, 0))
    group = list(permutations(range(4)))
    # On any free rank-three module the exterior square is the dual times det.
    J = s.Matrix([[0, 0, 1], [0, -1, 0], [1, 0, 0]])
    assert all(J*exterior(T) == T.det()*T.inv().T*J for T in map(rep, group))
    gram = s.eye(3)-s.ones(3)/4  # canonical quotient of the four equally weighted loops
    assert all(T.T*gram*T == gram for T in map(rep, group))
    x = s.Symbol('lambda')
    assert deck.charpoly().as_expr() == x**3+x**2+x+1
    assert exterior(deck).charpoly().as_expr() == s.expand((x-1)*(x**2+1))
    W = [exterior(rep(p)) for p in group]
    # Exact simultaneous intertwiner with the original v117 matrices, not just
    # a match of eigenvalues of the deck generator.
    M0 = exact_clock_dictionary()['M']
    U = s.diag(1, s.I, -s.I)
    assert exterior(M0) == exterior(-M0)
    C2 = J*exterior(M0)*J.inv()
    assert all(s.simplify(M.det()*C2.inv().T-M) == s.zeros(3) for M in [M0,-M0])
    lepton_sign_split = [s.simplify((s.eye(3)-M/2).det()) for M in [M0,-M0]]
    assert lepton_sign_split == [rat(7, 8), rat(9, 8)]
    cycle = s.Matrix(3, 3, lambda i, j: s.Symbol(f'z{i}{j}'))
    solutions = []
    for p, T in zip(group, W):
        if p == tuple(range(4)) or T**3 != s.eye(3):
            continue
        C, _ = s.linear_eq_to_matrix(list(cycle*exterior(deck)-U*cycle)
                                    + list(cycle*T-M0*cycle), list(cycle))
        kernel = C.nullspace()
        if kernel:
            assert len(kernel) == 1
            S = kernel[0].reshape(3, 3)
            assert s.simplify(S.det()) != 0
            solutions.append((p, S))
    assert len(solutions) == 4
    # Any intertwiner of the two original generators with themselves is scalar.
    C, _ = s.linear_eq_to_matrix(list(cycle*U-U*cycle)+list(cycle*M0-M0*cycle), list(cycle))
    assert C.nullspace() == [s.eye(3).reshape(9, 1)]
    # Orthogonal rank-one orbit in the two-form metric: full S4 removes all directions.
    b = s.Matrix(pi_left(K))
    G2dual = exterior(gram.inv())
    projector = b*(b.T*G2dual)/(b.T*G2dual*b)[0]
    cov_actions = [T.inv().T for T in W]
    twirl = sum((T*projector*T.inv() for T in cov_actions), s.zeros(3))/len(W)
    assert twirl == s.eye(3)/3

    rows = s.Matrix.vstack(ONE.T*K, A.T*K)
    full = rows.row_join(-rows*s.ones(3, 1))
    read = lambda inds: sum(abs(full[:, [i, j]].det()) for i in inds for j in inds if i < j)
    marks = [int(read(set(range(4))-{i})) for i in range(4)]
    assert marks == [19, 21, 11, 11] and read(range(4)) == 31
    # Under the literal H1 identification, one unobserved orthogonal channel
    # cannot generate K even up to overall scale: its largest pair would coincide.
    sing2 = (gram.inv()*K.T*gram*K).charpoly().as_expr()
    assert s.discriminant(sing2) != 0
    # Stronger: no identification with ANY positive D4 invariant metric in
    # the compiler lattice can make K a scaled one-channel compression.
    t = s.Symbol('t', positive=True)
    H = P.inv().T*s.diag(t, 1, 1)*P.inv()
    assert G.T*H*G == H and TA.T*H*TA == H
    pol = (H.inv()*K.T*H*K).charpoly()
    numerator = s.together(s.discriminant(pol.as_expr(), pol.gen)).as_numer_denom()[0]
    assert all(c > 0 for c in s.Poly(numerator, t).all_coeffs())
    anchor_projector = (s.eye(3)+TA)/2
    assert anchor_projector.rank() == 2 and anchor_projector*ONE == ONE and anchor_projector*A == A
    projected = [(anchor_projector*M*anchor_projector).adjugate().trace() for M in [Q,K,L]]
    assert projected == [rat(9, 2), -7, 0]
    return {'paper_ratios': list(map(str, actual)), 'left_Q_L_sums': [sum(pi_left(Q)), sum(pi_left(L))],
            'scalar_gcds': gcds, 'singlet_scalar_ud_roots': list(map(float, delta_roots)),
            'old_sector_products_e_u_d': old_products, 'active_word_stripped_det_ratio': expected_volume_ratio,
            'target_det_GdR_over_det_GuR_if_kinetic_normalization_only': target_right_gram,
            'delta_archived_first_order_pole_halfstep': delta_test, 'old_mobius_mu_tau': mobius,
            'new_scalar_ce_cmu': c6_lepton_e_mu,
            'four_left_mark_readouts': marks, 'full_left_readout': 31, 'full_twirl': str(twirl),
            'intertwiner_first_three_cycle_and_matrix': [solutions[0][0], str(solutions[0][1])],
            'intertwiner_number_of_cyclic_choices': len(solutions),
            'lepton_sheet_determinants_from_exterior_plus_orientation': lepton_sign_split,
            'K_metric_squared_singular_polynomial': str(sing2),
            'D4_K_metric_discriminant_numerator': str(numerator), 'anchor_projection_det_QKL': projected,
            'left_invariant_squared_norm': str((b.T*G2dual*b)[0]),
            'right_invariant_squared_norm': str((s.Matrix(pi_right(K)).T*exterior(gram)*s.Matrix(pi_right(K)))[0])}


if __name__ == '__main__':
    print(json.dumps({'exact': test(), 'ODE': [monodromy(.19, 2e-11), monodromy(.31, 5e-13)]}, default=str, indent=2))
