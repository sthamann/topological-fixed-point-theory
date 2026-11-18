"""Exact tests of what primitive spectral data determine about a response.

This tests the stated relative spectral formula, not a physical P1 insertion.
The 2x2 pair is a counterexample to an inference from finite channel values.
It is not substituted for the original geometric source.
"""
from pathlib import Path
import json
import mpmath as mp
import sympy as sp

HERE = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def exact_parity_variation():
    """A general real Hermitian 2x2 test of the even/odd variation identities."""
    a, b, z, u, v, w = sp.symbols("a b z u v w", real=True)
    r = sp.diag(1, -1)
    even = sp.diag(a, b)
    odd = sp.Matrix([[0, z], [z, 0]])
    he = sp.diag(u, v)
    ho = sp.Matrix([[0, w], [w, 0]])
    d = even + odd
    h = he + ho
    vd = -(h*d+d*h)/2
    va = (vd+r*vd*r)/2
    vb = (vd-r*vd*r)/2
    expected_a = -(he*even+even*he+ho*odd+odd*ho)/2
    expected_b = -(he*odd+odd*he+ho*even+even*ho)/2
    require(sp.simplify(va-expected_a) == sp.zeros(2), "even variation")
    require(sp.simplify(vb-expected_b) == sp.zeros(2), "odd variation")
    return {"split_variation_exact": True,
            "reference_variation": "-({h_even,A}+{h_odd,B})/2",
            "relative_variation": "-({h_even,B}+{h_odd,A})/2"}


def same_value_opposite_force():
    mp.mp.dps = 90
    rational_q = [sp.Rational(1, 2), sp.Rational(2, 3),
                  sp.Rational(3, 4), sp.Rational(4, 5)]
    moment_matrix = sp.Matrix([[1]*4, rational_q, [q**4 for q in rational_q]])
    null = moment_matrix.nullspace()[0]
    null /= max(abs(x) for x in null)
    require(moment_matrix*null == sp.zeros(3, 1), "exact channel matching")
    conv = lambda x: mp.mpf(str(sp.numer(x)))/mp.mpf(str(sp.denom(x)))
    q = [conv(x) for x in rational_q]
    t = [-mp.log(x) for x in q]
    c = [tj*(qj-4*qj**4) for tj, qj in zip(t, q)]
    # Exact nonzero certificate using independence of logarithms of primes.
    log2,log3,log5 = sp.symbols("log2 log3 log5")
    log_q = [-log2,log2-log3,log3-2*log2,2*log2-log5]
    exact_force_difference = sp.expand(sum(-lq*(qj-4*qj**4)*vj
                                         for lq,qj,vj in zip(log_q,rational_q,null)))
    require(exact_force_difference != 0,"nonzero rational prime-log combination")
    # A separate outward interval computation certifies the sign.
    mp.iv.dps = 45
    qi = [mp.iv.mpf(int(sp.numer(x)))/int(sp.denom(x)) for x in rational_q]
    ni = [mp.iv.mpf(int(sp.numer(x)))/int(sp.denom(x)) for x in null]
    difference_interval = sum(-mp.iv.log(x)*(x-4*x**4)*v for x,v in zip(qi,ni))
    require(bool(difference_interval > 0),"certified positive force difference")
    # Choose a strictly positive center measure with exactly zero scale force.
    weights = [-sum(c[1:])/c[0], mp.mpf(1), mp.mpf(1), mp.mpf(1)]
    normalization = sum(weights)
    weights = [x/normalization for x in weights]
    delta = min(weights)/2
    null_numeric = [conv(x) for x in null]
    plus = [x+delta*v for x, v in zip(weights, null_numeric)]
    minus = [x-delta*v for x, v in zip(weights, null_numeric)]
    require(min(plus+minus) > 0, "strictly positive profile measures")

    def profile(measure, eigenvalue):
        return sum(wj*mp.exp(-tj*eigenvalue**2) for wj,tj in zip(measure,t))

    def action(measure, scale):
        ratio = mp.exp(-scale)
        return 2*(profile(measure, ratio)-profile(measure, 2*ratio))

    # A=2I, B=sigma_x, R=diag(1,-1), D=A+B. Scale BOTH pieces.
    r = sp.diag(1,-1)
    a = 2*sp.eye(2)
    b = sp.Matrix([[0,1],[1,0]])
    d = a+b
    require((d+r*d*r)/2 == a and (d-r*d*r)/2 == b, "original split")
    derivative_plus = 4*sum(wj*cj for wj,cj in zip(plus,c))
    derivative_minus = 4*sum(wj*cj for wj,cj in zip(minus,c))
    exact_equal_error = max(abs(profile(plus,x)-profile(minus,x)) for x in (0,1,2))
    require(exact_equal_error < mp.mpf("1e-85"), "same operator channel values")
    require(derivative_plus > 0 and derivative_minus < 0, "opposite force")
    require(abs(derivative_plus+derivative_minus) < mp.mpf("1e-85"), "opposite centered slopes")
    derivative_error = max(abs(mp.diff(lambda s: action(ws,s),0)-slope)
                           for ws,slope in ((plus,derivative_plus),(minus,derivative_minus)))
    require(derivative_error < mp.mpf("1e-85"), "independent differentiated action")
    # The eta phase is held fixed in this algebraic spectral test. No eigenvalue
    # of A,B,D crosses zero on the scale path. No physical boundary eta is inferred.
    return {"status": "exact_scoped_nonidentifiability",
            "pair": {"A":"2 I", "B":"sigma_x", "R":"diag(1,-1)",
                     "D":"A+B", "path":"D(s)=exp(-s)D; chi=1"},
            "q": [str(x) for x in rational_q],
            "heat_times": [float(x) for x in t],
            "exact_null_vector": [str(x) for x in null],
            "exact_prime_log_difference":str(exact_force_difference),
            "positive_difference_interval":str(difference_interval),
            "null_vector_annihilates_rows": ["1","q","q^4"],
            "weights_plus": [float(x) for x in plus],
            "weights_minus": [float(x) for x in minus],
            "matched_profile_values": {str(x):float(profile(plus,x)) for x in (0,1,2)},
            "same_action_at_s0": float(action(plus,0)),
            "derivative_plus": float(derivative_plus),
            "derivative_minus": float(derivative_minus),
            "profile_equality_error_90_digits": float(exact_equal_error),
            "independent_derivative_error_90_digits": float(derivative_error),
            "not_implied": ["either measure is the original TFPT transfer profile",
                            "finite 2x2 example is the P1 geometry",
                            "a specified full transfer channel permits both profiles"],
            "curves": [{"s":j/100,
                        "S_plus":float(action(plus,mp.mpf(j)/100)),
                        "S_minus":float(action(minus,mp.mpf(j)/100))}
                       for j in range(-25,26)]}


def coupled_force_check():
    """Independently differentiate the full matrix expression with moving chi."""
    mp.mp.dps = 90
    r = mp.matrix([[1,0],[0,-1]])
    d = mp.matrix([[2,mp.mpf(1)/3],[mp.mpf(1)/3,3]])
    h = mp.matrix([[mp.mpf(1)/5,mp.mpf(2)/7],
                   [mp.mpf(2)/7,-mp.mpf(1)/4]])
    a,b = (d+r*d*r)/2,(d-r*d*r)/2
    he,ho = (h+r*h*r)/2,(h-r*h*r)/2
    t,chi,kappa = mp.mpf(3)/5,mp.mpf(7)/6,mp.mpf(1)/9
    tr = lambda m:sum(m[j,j] for j in range(m.rows))
    ea,eb = mp.expm(-t*a*a/chi**2),mp.expm(-t*b*b/chi**2)
    even_force = 2*t/chi**2*tr(he*(b*b*eb-a*a*ea))
    odd_force = t/chi**2*tr(ho*(a*b*eb+b*eb*a-b*a*ea-a*ea*b))
    scale_force = 2*kappa*t/chi**2*tr(b*b*eb-a*a*ea)

    def action(s):
        factor = mp.expm(-s*h/2)
        ds = factor*d*factor
        ac,bc = (ds+r*ds*r)/2,(ds-r*ds*r)/2
        cs = chi*mp.exp(kappa*s)
        return tr(mp.expm(-t*bc*bc/cs**2)-mp.expm(-t*ac*ac/cs**2))

    derivative = mp.diff(action,0)
    error = abs(derivative-even_force-odd_force-scale_force)
    require(error < mp.mpf("1e-85"),"coupled force with odd profile and chi variation")
    return {"even_force":float(even_force),"odd_force":float(odd_force),
            "chi_force":float(scale_force),"full_derivative":float(derivative),
            "error_90_digits":float(error),
            "scope":"finite independent check of the general trace identity; eta held fixed"}


def universal_scale_stationarity():
    """Proof certificate: Laplace exponentials of distinct positive values separate."""
    # For finite positive squared spectra, universal stationarity for every
    # atomic positive heat measure implies sum multiplicity*x*exp(-t*x)=0.
    # Linear independence forces equal multiplicities at each x>0. The
    # zero-eigenvalue counts are invisible to this derivative and stay constant.
    x1,x2,x3,t = sp.symbols("x1 x2 x3 t", positive=True)
    wronskian = sp.wronskian([sp.exp(-t*x) for x in (x1,x2,x3)],t)
    expected = -(x2-x1)*(x3-x1)*(x3-x2)*sp.exp(-t*(x1+x2+x3))
    require(sp.simplify(wronskian-expected)==0,"independent heat responses")
    return {"finite_spectrum_theorem":
            "For a fixed-chi common scaling, stationarity for every positive atomic heat profile forces identical nonzero squared spectral multisets; the action is then constant on that entire scale orbit.",
            "zero_eigenvalues": "their multiplicities are not fixed by scale stationarity",
            "scope": "finite selfadjoint pair; fixed chi; full common scale path admitted; eta phase fixed; profiles not otherwise fixed by original channel",
            "wronskian_three_distinct_values_exact": True}


def retained_scalar_channel():
    q,alpha,delta = sp.symbols("q alpha delta",positive=True)
    b = sp.Rational(5,4)
    log_t = -b*sp.log(1-q)
    derivative_q = sp.diff(log_t,q)
    require(sp.simplify(derivative_q-b/(1-q))==0,"fixed channel derivative")
    original_q = delta*sp.exp(-2*alpha)
    derivative_alpha = derivative_q.subs(q,original_q)*sp.diff(original_q,alpha)
    require(sp.simplify(derivative_alpha+2*b*original_q/(1-original_q))==0,
            "original conditional alpha derivative")
    return {"given_channel":"T(q)=(1-q)^(-5/4), 0<q<1",
            "d_logT_dq":"5/[4(1-q)]",
            "given_old_weight":"q(alpha)=delta_top exp(-2 alpha)",
            "d_logT_dalpha":"-(5/2)q/(1-q)",
            "not_missing":"another arbitrary spectral profile for this fixed scalar channel",
            "still_required":"actual geometric insertion map to delta q and all retained residual/reference/measure terms"}


def main():
    result = {"parity_variation":exact_parity_variation(),
              "same_value_opposite_force":same_value_opposite_force(),
              "coupled_force_check":coupled_force_check(),
              "retained_scalar_channel":retained_scalar_channel(),
              "universal_scale_stationarity":universal_scale_stationarity()}
    path = HERE/"runs"/"primitive_selection.json"
    path.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({"result":str(path),"status":"exact_scoped",
                      "same_action":result["same_value_opposite_force"]["same_action_at_s0"],
                      "slopes":[result["same_value_opposite_force"]["derivative_plus"],
                                result["same_value_opposite_force"]["derivative_minus"]]},indent=2))


if __name__ == "__main__":
    main()
