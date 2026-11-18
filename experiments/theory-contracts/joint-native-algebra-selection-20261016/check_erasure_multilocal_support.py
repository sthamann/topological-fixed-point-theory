"""Literal native erasure20 and joint25 highest-CAR support under the original fold.

Preserves partial_a P_i / sqrt(N_i), with P_i the raw original quartic,
original M/O and both quarters.
No new physical source choice, operation menu, or source-energy assertion.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import sys

BASE = Path(__file__).resolve().parent
import check_native_quartic_multilocal_support as q

s, m, old = q.s, q.m, q.old


def hashes():
    return {n: hashlib.sha256((BASE/n).read_bytes()).hexdigest() for n in (
        'check_native_quartic_multilocal_support.py',
        'NATIVE-QUARTIC-MULTILOCAL-SUPPORT.txt',
        'native_quartic_multilocal_support_results.json')}


def main(output):
    before = hashes()
    m.need(hashlib.sha256(q.M_PATH.read_bytes()).hexdigest() == q.M_SHA,
           'Original M source retained')
    ld = m._lattice_data()
    T = s.Matrix.hstack(*map(s.Matrix,m._TARGET_SIMPLE))
    O = old.cl(T*ld['basis'].inv()/s.sqrt(2))
    m.need(hashlib.sha256(json.dumps(q.rows(O),separators=(',',':')).encode()).hexdigest()
           == q.O_ROWS_SHA, 'Original O retained')
    M = s.Matrix(json.loads(q.M_PATH.read_text())['explicit_Weyl_matrix'])
    Uplus = s.zeros(8,4)
    for a in range(4):
        Uplus[2*a,a] = 1/s.sqrt(2)
        Uplus[2*a+1,a] = -s.I/s.sqrt(2)
    V = old.cl(M.T*O*Uplus)
    polys, norms = m.code()
    labels = [(i,a) for i in range(5) for a in range(4)]
    erasure = [m.deriv(polys[i], [m.Q.one if a==b else m.Q.zero for b in range(4)])
               for i,a in labels]
    original_gram = s.Matrix(20,20,lambda j,k:old.Q.to_sympy(m.inner(erasure[j],erasure[k])))
    weights = [norms[i] for i,a in labels]
    m.need(original_gram == s.diag(*weights), 'Actual twenty unit erasure fields are orthonormal')
    subsets = list(itertools.combinations(range(8),3))
    C = s.Matrix([
        [s.simplify(sum(old.Q.to_sympy(c)*sum(
            s.prod(V[k,a] for k,a in zip(K,perm))
            for perm in itertools.permutations(word)) for word,c in f.items()))
         for f in erasure] for K in subsets])
    G = old.cl(C.conjugate().T*C)
    m.need(C.shape == (56,20) and C.rank() == 20, 'Literal cubic principal map injective')
    N = s.diag(*[1/s.sqrt(n) for n in weights])
    Gn = old.cl(N*G*N)
    D = s.Matrix.vstack(*[C[j:j+1,:] for j,K in enumerate(subsets) for _ in K])
    Ga = old.cl(D.conjugate().T*D/8)
    m.need(Ga == s.Rational(3,8)*G and Ga.rank() == 20,
           'Every original sheet has injective degree-five CAR commutator')
    h = s.symbols('h0:4')
    hpoly = [sum(old.Q.to_sympy(c)*s.prod(h[a] for a in word) for word,c in f.items())
             for f in erasure]
    quarters = {}
    for alpha in (s.Rational(1,4),s.Rational(3,4)):
        delta = s.Matrix([(k+alpha)/8-s.Rational(1,2) for k in range(8)])
        b = old.cl(V.T*delta)
        for f in hpoly:
            p = s.Poly(s.expand(f.subs(dict(zip(h,[h[a]+b[a]for a in range(4)])),
                                     simultaneous=True)),*h)
            top = sum(c*s.prod(h[a]**mon[a]for a in range(4))
                      for mon,c in p.terms()if sum(mon)==3)
            m.need(s.expand(top-f)==0,'Original b shift preserves cubic principal term')
        quarters[str(alpha)] = {'b':[str(x)for x in b], 'cubic_shift_unchanged':True,
                                'all_sheet_ranks_from_uniform_gram_identity':[20]*8,
                                'sheet_gram_factor':'3/8'}
    _, C4 = q.principal_coefficients(polys,V)
    m.need(C4.rank()==5, 'Retain reviewed quartic principal rank')
    after = hashes()
    m.need(before==after,'All reviewed source/result files untouched')
    result = {
        'status':'PASS_SCOPED',
        'source_inputs':{'M_file_sha256':q.M_SHA,'O_canonical_rows_sha256':q.O_ROWS_SHA,
                         'V':q.rows(V)},
        'column_order':[[i,a]for i,a in labels],
        'field_definition':'g_(i,a)=partial_a P_i / sqrt(N_i), P_i raw original quartic',
        'original_norms':weights, 'original_gram':q.rows(original_gram),
        'subset_order':[list(K)for K in subsets],
        'principal_coefficients_raw':q.rows(C), 'principal_gram_raw':q.rows(G),
        'principal_rank':C.rank(), 'principal_nullity':20-C.rank(),
        'principal_gram_raw_determinant':str(G.det()),
        'principal_gram_unit':q.rows(Gn),
        'principal_gram_unit_determinant':str(Gn.det()),
        'sheet_commutator_gram_raw':q.rows(Ga),
        'sheet_commutator_rank':Ga.rank(), 'quarters':quarters,
        'joint25_graded_rank':C4.rank()+C.rank(),
        'joint25_reason':'Any nonzero quartic coefficient survives CAR degree8 / commutator7; when it vanishes, any nonzero erasure coefficient survives degree6 / commutator5. No cancellation across these highest different filtration degrees.',
        'scope':'Given original 8-fold M/O current dictionary only; no P1 source/reference/preparation selection or all-operations no-go.',
        'protected_files_before':before, 'protected_files_after':after}
    output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:result[k]for k in ('status','principal_rank','principal_nullity',
                                         'principal_gram_raw_determinant',
                                         'principal_gram_unit_determinant',
                                         'sheet_commutator_rank','joint25_graded_rank')}))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--output',type=Path,default=BASE/'erasure_multilocal_support_results.json')
    main(p.parse_args().output)
