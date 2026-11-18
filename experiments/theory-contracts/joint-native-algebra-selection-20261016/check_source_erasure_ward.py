"""Exact source-erasure Ward identity and uncompressed return leakage.

Uses the original quartic code, actual erasure derivative and previously
derived full lattice vertex action. No code projection is applied to an
intermediate state. All checks remain active under Python -OO.
"""
from pathlib import Path
import json
import sympy as s
from check_mixed_return import Q, need, conj, add, mul, deriv, linear, inner, code, response
import check_mixed_return as original


def annihilate(p, a):
    """bar(u_a)(1) on charged Fock output; other modes stay untouched."""
    out = {}
    for mon, value in p.items():
        for j, k in enumerate(mon):
            if k == a:
                n = mon[:j] + mon[j+1:]
                out[n] = out.get(n, Q.zero) + value
    return {m: c for m, c in out.items() if c}


def contact_modes(f, b, c):
    """E_beta,1 f and W_beta,2 f from original Laurent coefficients."""
    B1, B2 = linear(b, 1), linear(b, 2)
    V1, V2 = linear(c, 1), linear(c, 2)
    E2 = add((Q(1)/2, mul(B1, B1)), (Q(1)/2, B2))
    F = [f]
    for k in range(1, 5):
        F.append(add((-Q.one/k, deriv(F[-1], b))))
    G = [deriv(f, c)]
    for k in range(1, 4):
        G.append(add((-Q.one/k, deriv(G[-1], b))))
    e1 = add((Q.one, F[2]), (Q.one, mul(F[3], B1)), (Q.one, mul(F[4], E2)))
    w2 = add((Q.one, mul(V1, add((Q.one, F[3]), (Q.one, mul(F[4], B1))))),
             (Q.one, mul(V2, F[4])), (Q.one, G[1]),
             (Q.one, mul(G[2], B1)), (Q.one, mul(G[3], E2)))
    return e1, w2


def main():
    polys, norms = code()
    b = [Q.zero, Q.zero, Q(0, -1), Q.zero]
    c = [Q.zero, Q.zero, Q.zero, Q(4, -4)]
    source = [response(f, b, c) for f in polys]
    checks, nonzero_contacts = 0, 0
    for j, f in enumerate(polys):
        e1, w2 = contact_modes(f, b, c)
        need(all(sum(1+k//8 for k in m) == 2 for p in (e1, w2) for m in p),
             'Contact outputs have oscillator grade2 / total grade3 and charge beta')
        for a in range(4):
            da = [Q.one if k == a else Q.zero for k in range(4)]
            df = deriv(f, da)
            forward = response(df, b, c) if df else {}
            actual = add((Q.one, annihilate(source[j], a)), (-Q.one, forward))
            contact = add((conj(c[a]), e1), (conj(b[a]), w2))
            need(actual == contact, f'Exact source-erasure Ward identity {j},{a}')
            checks += 1
            nonzero_contacts += bool(contact)
    # h=x0^4/sqrt24 is killed by the ACTUAL source operation.
    h = {(0, 0, 0, 0): Q.one}
    need(not response(h, b, c), 'Original source annihilates h exactly')
    q = add((Q.one, h), (-Q(1)/4, polys[0]))
    need(inner(q, q) == Q(18), 'q_raw norm squared18')
    for f in polys:
        need(inner(f, q) == Q.zero, 'q is orthogonal to all5 original code states')
    Sq = response(q, b, c)
    need(Sq == add((-Q(1)/4, source[0])), 'Exact image of code-orthogonal witness')
    raw_gram00 = Q.to_sympy(inner(source[0], source[0]))
    need(raw_gram00 == 6016, 'Original full-source response norm6016')
    overlap = s.simplify(Q.to_sympy(inner(Sq, source[0])) / s.sqrt(18*norms[0]))
    need(s.simplify(overlap + 188/(3*s.sqrt(3))) == 0, 'Nonzero UNCOMPRESSED return leakage')
    leak_bound = s.simplify(s.conjugate(overlap)*overlap)
    need(leak_bound == s.Rational(35344, 27), 'Exact norm-squared leakage lower bound')
    # Independent elementary amplitude calculation, without q source evaluation.
    K00 = s.Rational(raw_gram00, norms[0])
    need(s.simplify(-K00/s.sqrt(3) - overlap) == 0, 'Independent annihilated-h argument')
    need(not annihilate(source[0], 0), 'Erasure a0 of actual forward response is zero')
    need(not response(deriv(polys[0], [Q.one, Q.zero, Q.zero, Q.zero]), b, c),
         'Actual field after erasure a0 is zero')
    # Original fixed root-input-basis sum, with NO new weights or selected root.
    lattice = original._lattice_data()
    T = s.Matrix.hstack(*map(s.Matrix, original._TARGET_SIMPLE))
    O = original.clean(T*lattice['basis'].inv()/s.sqrt(2))
    native = original._native_neutral_branch()
    source_roots = original._actual_current_source()['chevalley'].roots
    ri = original.clean(s.sqrt(2)*O.T/2)
    sums = []
    for omitted in range(4):
        subset = [r for r in native['spinor'] if tuple(r[5:]) == tuple(native['family'][omitted])]
        A16 = sum((s.Matrix(r)*s.Matrix(r).T/4 for r in subset), s.zeros(8))
        AN = original.clean(O.T*(4*s.eye(8)-A16)*O)
        aggregate = Q.zero
        for source_root in source_roots:
            r = ri*s.Matrix(source_root)
            need(all(t.is_Integer for t in r), 'Original root map has integral native roots')
            vv = -AN*r + ((r.T*AN*r)[0]/4)*r
            rb = [Q.from_sympy((r[2*j]-s.I*r[2*j+1])/2) for j in range(4)]
            rc = [Q.from_sympy((vv[2*j]-s.I*vv[2*j+1])/2) for j in range(4)]
            if omitted == 0 and tuple(source_root) == (0,0,0,0,0,0,-2,-2):
                need(rb == b and rc == c, 'Witness transported from actual original source')
            aggregate += inner(response(q, rb, rc), response(polys[0], rb, rc))
        normalized = s.simplify(Q.to_sympy(aggregate)/s.sqrt(18*norms[0]))
        sums.append({'omitted_family': omitted, 'raw_overlap': str(Q.to_sympy(aggregate)),
                     'normalized_off_code_amplitude': str(normalized),
                     'leakage_norm_squared_lower_bound': str(s.simplify(normalized*s.conjugate(normalized)))})
        print('original family', omitted, 'aggregate off-code amplitude', normalized, flush=True)
    need(sums[0]['raw_overlap'] == '-50208', 'Original aggregate off-code witness')
    need(all(x['raw_overlap'] == '-50208' and x['leakage_norm_squared_lower_bound'] == '4376464/3'
             for x in sums), 'All4 original family omissions have the same nonzero leakage witness')
    out = {
        'status': 'PASS_EXACT_SOURCE_ERASURE_WARD_AND_UNCOMPRESSED_CODE_LEAKAGE',
        'contract_id': 'UR.NATIVE.JOINT_OPERATOR.01',
        'source_root_doubled': [0,0,0,0,0,0,-2,-2],
        'b': [str(Q.to_sympy(x)) for x in b],
        'c': [str(Q.to_sympy(x)) for x in c],
        'mode_convention': 'Y(A,z)=sum_n A_n z^(-n-1); o(W)=W_1 for wt(W)=2',
        'source_erasure_ward': '[bar(u_a)(1),W_1]=conj(c_a) E_beta,1 + conj(b_a) W_beta,2',
        'exact_code_ward_identities': checks,
        'nonzero_contact_instances': nonzero_contacts,
        'raw_code0_forward_norm_squared': str(raw_gram00),
        'K00': str(K00),
        'h': 'x0^4/sqrt24; <h,f0>=1/2; W_1 h=0',
        'normalized_code_orthogonal_q': '(x0^4-(x0^4+x1^4+x2^4+x3^4)/4)/sqrt18',
        'uncompressed_return_off_code_amplitude': str(overlap),
        'uncompressed_return_leakage_norm_squared_lower_bound': str(leak_bound),
        'projection_artifact': '<h,L f0>=0 but <h,Pi L Pi f0>=94/3',
        'a0_contact_coefficient': '0 for both E_beta,1 and W_beta,2',
        'original_root_input_sum': 'sum_beta o(W_beta)^dagger o(W_beta); no inter-root cross terms or reweighting. This is the original input-basis response sum, not a derived physical source selection.',
        'original_family_aggregate_leakage': sums,
        'scope': 'Given original E8_1 source, native positive Cartan embedding and supplied R_neut branch. This specific full positive adjoint return does not preserve the code. Endpoint-compressed K5 is not the full code operation. No physical instrument, source history, P1 field selection or3+2 selection is derived. No exclusion of all source histories.'
    }
    Path(__file__).with_name('source_erasure_ward_results.json').write_text(json.dumps(out, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(out, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
