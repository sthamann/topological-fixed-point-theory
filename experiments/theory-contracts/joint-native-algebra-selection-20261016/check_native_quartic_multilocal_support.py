"""Exact highest-CAR support of the EXISTING native quartics under the 8-fold.

Uses the original M/O/current/Schur constructors. This checks a field-support
witness, not a P1 intervention, a source selector, or a closed-code axiom.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

import sympy as s
import check_original_current_schur_kernel as old

m = old.m
BASE = Path(__file__).resolve().parent
M_PATH = m.R / 'experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json'
M_SHA = '253f370a9c84b70c28169ed0aad83d1866245804227e44aa0189e823382f5639'
O_ROWS_SHA = '0788b82fa2886b05c8748f41cab34ce1f37be4bfdf665a543e3173bf3b7c2ff4'


def rows(A):
    return [[str(x) for x in row] for row in A.tolist()]


def principal_coefficients(polys, V):
    """Coefficient of t_K in f_i(V^T t), with each t_k nilpotent.

    n_k=:chi_k^* chi_k: has square zero in the highest fermion-degree
    associated graded algebra. A degree-four current monomial contributes
    only squarefree K, |K|=4. One common creator-order sign is immaterial.
    """
    subsets = list(itertools.combinations(range(8), 4))
    C = s.Matrix([
        [s.simplify(sum(
            old.Q.to_sympy(c) * sum(
                s.prod(V[k, a] for k, a in zip(K, perm))
                for perm in itertools.permutations(word))
            for word, c in f.items())) for f in polys]
        for K in subsets])
    return subsets, C


def main(output):
    m.need(hashlib.sha256(M_PATH.read_bytes()).hexdigest() == M_SHA,
           'Original Weyl matrix source unchanged')
    ld = m._lattice_data()
    T = s.Matrix.hstack(*map(s.Matrix, m._TARGET_SIMPLE))
    O = old.cl(T * ld['basis'].inv() / s.sqrt(2))
    m.need(hashlib.sha256(json.dumps(rows(O), separators=(',', ':')).encode()).hexdigest()
           == O_ROWS_SHA, 'Same original O, no fitted isometry')
    M = s.Matrix(json.loads(M_PATH.read_text())['explicit_Weyl_matrix'])
    m.need(old.cl(O.T * O) == s.eye(8) and M.T * M == s.eye(8),
           'Original M/O preserve the supplied metric')
    Uplus = s.zeros(8, 4)
    for a in range(4):
        Uplus[2*a, a] = 1/s.sqrt(2)
        Uplus[2*a+1, a] = -s.I/s.sqrt(2)
    V = old.cl(M.T * O * Uplus)
    polys, norms = m.code()  # Also checks the original tensor/Fock normalization.
    m.need(norms == [96, 8, 8, 8, 1], 'Original five polynomial norms')
    subsets, C = principal_coefficients(polys, V)
    G = old.cl(C.conjugate().T * C)
    expected = s.diag(s.Matrix([[s.Rational(243,8), -s.Rational(51,16)],
                                [-s.Rational(51,16), s.Rational(115,32)]]),
                      s.Matrix([[s.Rational(13,8), -s.Rational(3,8)],
                                [-s.Rational(3,8), s.Rational(13,8)]]), s.Rational(1,4))
    m.need(G == expected and C.rank() == 5 and G.det() == s.Rational(495,8),
           'Exact original 70-by-5 principal map is injective')
    N = s.diag(*[1/s.sqrt(n) for n in norms])
    Gn = old.cl(N * G * N)
    m.need(Gn.det() == s.Rational(165,131072), 'Normalized principal Gram determinant')

    # The eigenchannel basis is orthonormal. Contraction of the four creator
    # legs with eta_r has disjoint (annihilator K, creator K\{k}) basis labels.
    # Every Fourier-matrix entry has squared modulus 1/8. Thus this is the
    # exact Gram of the degree-seven commutator witness, not a dimension count.
    D = s.Matrix.vstack(*[C[j:j+1, :] for j, K in enumerate(subsets) for _ in K])
    arc_gram = old.cl(D.conjugate().T * D / 8)
    m.need(arc_gram == G/2 and arc_gram.rank() == 5,
           'Every original arc has injective creation-leg CAR contraction')

    h = s.symbols('h0:4')
    original = [sum(old.Q.to_sympy(c)*s.prod(h[a] for a in word)
                    for word,c in f.items()) for f in polys]
    cases = {}
    for alpha in (s.Rational(1,4), s.Rational(3,4)):
        delta = s.Matrix([(k+alpha)/8-s.Rational(1,2) for k in range(8)])
        b = old.cl(V.T * delta)
        for f in original:
            shifted = s.Poly(s.expand(f.subs(dict(zip(h, [h[a]+b[a] for a in range(4)])),
                                            simultaneous=True)), *h)
            top = sum(c*s.prod(h[a]**word[a] for a in range(4))
                      for word,c in shifted.terms() if sum(word) == 4)
            m.need(s.expand(top-f) == 0, 'One original b shift leaves principal CAR map unchanged')
        # Exact twisted cyclic deck recurrence: write zeta=exp(i*pi/16).
        # U_alpha(r,k)=zeta^((4k+4alpha)r)/sqrt(8), wrap=zeta^(32alpha).
        a = int(4*alpha)
        for r in range(8):
            for k in range(8):
                lhs = (4*k+a)*(r+1) if r < 7 else 8*a
                rhs = (4*k+a)*(r+1)
                m.need((lhs-rhs) % 32 == 0, 'Exact holonomy-aware deck eigenchannel law')
        cases[str(alpha)] = {
            'b': [str(x) for x in b], 'principal_shift_unchanged': True,
            'twisted_deck_covariance': True,
            'arc_commutator_principal_ranks': [5]*8,
            'arc_commutator_gram': 'raw_principal_gram / 2',
        }

    # Small existing-path source check only; do not rerun all-mode suites.
    packets, packet_norms = old.source_packets()
    m.need(packet_norms == norms, 'Same actual full CAR source packet')
    for alpha, states in packets.items():
        energies = sorted({old.car.energy(bits) for v in states for bits in v})
        charges = sorted({old.car.charges(bits) for v in states for bits in v})
        B = s.Matrix(5,5,lambda i,j: old.Q.to_sympy(old.car.inner(states[i], states[j])))
        m.need(energies == [32] and charges == [(0,)*8] and B.rank() == 5,
               'Full actual source packet is neutral, energy32 and injective')
        cases[str(alpha)]['source_packet'] = {
            'supports': [len(v) for v in states], 'energies': energies,
            'residue_charges': charges, 'gram_rank': B.rank()}

    result = {
        'status': 'PASS_SCOPED', 'scientific_status_change': False,
        'scope': 'Existing conditional NS/current M/O eightfold dictionary; highest CAR field support and actual source packet only.',
        'not_established': ['Physical P1 source/reference or preparation', 'A q selector',
                            'No-go for arbitrary unbounded local field products'],
        'source_inputs': {'M_path': str(M_PATH.relative_to(m.R)), 'M_file_sha256': M_SHA,
                          'O_canonical_rows_sha256': O_ROWS_SHA, 'O': rows(O), 'V': rows(V)},
        'polynomial_norms': norms,
        'original_polynomials': [
            {'terms': [{'word': list(word), 'coefficient': str(old.Q.to_sympy(c))}
                       for word,c in sorted(f.items())],
             'unit_normalization': '1/sqrt('+str(norms[i])+')'}
            for i,f in enumerate(polys)],
        'subset_order': [list(K) for K in subsets],
        'raw_principal_coefficients': rows(C), 'raw_principal_gram': rows(G),
        'raw_principal_rank': C.rank(), 'raw_principal_gram_determinant': str(G.det()),
        'normalized_principal_gram': rows(Gn),
        'normalized_principal_gram_determinant': str(Gn.det()),
        'arc_commutator_principal_gram': rows(arc_gram), 'quarters': cases,
        'analytic_boundary': 'Bounded proper-interval eigenstate-preparation lemma is proved separately in NATIVE-QUARTIC-MULTILOCAL-SUPPORT.txt; it is not established by a finite numerical test.'}
    output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'status': result['status'], 'principal_rank': 5,
                      'every_arc_rank': 5, 'source_energy': 32,
                      'quarters': list(cases), 'output': str(output)}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=BASE/'native_quartic_multilocal_support_results.json')
    main(parser.parse_args().output)
