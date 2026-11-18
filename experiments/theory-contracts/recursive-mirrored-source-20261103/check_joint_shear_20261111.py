"""Exact local test of the two response directions at fixed flat sources."""
import json
import sympy as s


def run():
    x, y, w, z = s.symbols("x y w z", real=True)
    eps, a, h, kg, zh, lg, lh = s.symbols("eps a h K_gamma Z_H lambda_gamma lambda_H", real=True)
    p, m = s.symbols("p m", positive=True)
    coords = (x, y, w, z)
    f = s.cos(2*m*z)
    g = s.diag(s.exp(-2*eps*f), s.exp(2*eps*f), 1, 1)
    ginv = g.inv()
    assert s.simplify(g.det()) == 1
    n = len(coords)
    gamma = {}
    for i in range(n):
        for j in range(n):
            for k in range(n):
                gamma[i,j,k] = s.simplify(sum(ginv[i,l] * (
                    s.diff(g[l,j], coords[k])+s.diff(g[l,k], coords[j])-
                    s.diff(g[j,k], coords[l])) for l in range(n))/2)
    ricci = s.Matrix(n,n,lambda i,j: s.simplify(sum(
        s.diff(gamma[k,i,j],coords[k])-s.diff(gamma[k,i,k],coords[j])
        +sum(gamma[k,k,l]*gamma[l,i,j]-gamma[k,j,l]*gamma[l,i,k] for l in range(n))
        for k in range(n))))
    ricci_jet = ricci.diff(eps).subs(eps,0).applyfunc(s.simplify)
    assert ricci.subs(eps,0) == s.zeros(4)
    assert ricci_jet == s.diag(s.diff(f,z,2),-s.diff(f,z,2),0,0)
    assert s.simplify(s.trace(ginv@ricci).diff(eps).subs(eps,0)) == 0
    field = s.cos(p*x)*s.cos(m*z)
    photon = [0, 0, a*field, 0]
    F = s.Matrix(n,n,lambda i,j: s.diff(photon[j],coords[i])-s.diff(photon[i],coords[j]))
    H = s.Matrix([s.diff(h*field,v) for v in coords])
    Fup = ginv@F@ginv
    D = ginv@H
    base = kg*sum(F[i,j]*Fup[i,j] for i in range(n) for j in range(n))/4 + zh*(H.T@D)[0]
    curvature = lg*s.trace(ricci@Fup@g@Fup.T) + lh*(D.T@ricci@D)[0]
    assert curvature.subs(eps,0) == 0
    assert s.simplify(curvature.subs({a:0,h:0})) == 0
    jets = [(s.diff(base+curvature,eps,v,v).subs(eps,0)).trigsimp() for v in (a,h)]
    common = p**2*s.sin(p*x)**2*s.cos(m*z)**2*f
    target = [2*kg-8*lg*m**2,4*zh-8*lh*m**2]
    assert all(s.simplify(j-common*t) == 0 for j,t in zip(jets,target))
    # Nonzero momentum-conserving mean: no spuriously nonzero integral of a
    # single nonconstant graviton with two zero longitudinal field momenta.
    avg = s.integrate(s.integrate(common,(x,0,2*s.pi/p)),(z,0,2*s.pi/m))*p*m/(4*s.pi**2)
    assert avg == p**2/8
    matrix = s.Matrix(target).jacobian([lg,lh])
    assert matrix.det() == 64*m**4
    assert matrix.subs(m,0) == s.zeros(2)
    m1, m2 = s.symbols("m1 m2", positive=True)
    intercepts = [s.factor((m2**2*t.subs(m,m1)-m1**2*t.subs(m,m2))/(m2**2-m1**2))
                  for t in target]
    assert intercepts == [2*kg, 4*zh]
    # Gauge-shift check on the exact field strength.
    gauge = s.Function("chi")(*coords)
    shifted = [photon[i]+s.diff(gauge,coords[i]) for i in range(n)]
    assert s.Matrix(n,n,lambda i,j: s.diff(shifted[j],coords[i])-s.diff(shifted[i],coords[j])) == F
    print(json.dumps({"normalized_shared_shear_response":[str(t) for t in target],
                      "curvature_coefficient_determinant":str(matrix.det()),
                      "momentum_matched_average":str(avg),
                      "soft_limit_fixed_coefficients": [str(t.subs(m,0)) for t in target],
                      "two_momentum_intercepts":[str(t) for t in intercepts],
                      "coefficient_response_abs_bound_for_abs_lambda_at_most_L":"2 L q^2, q=2m"},
                     indent=2))


if __name__ == "__main__":
    run()
