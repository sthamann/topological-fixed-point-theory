"""Recheck the manual radial derivatives and the quarter reference product."""

import mpmath as mp
import sympy as s

x, y, xx, xy, yy, v, lam, z, w = s.symbols('x y Gxx Gxy Gyy v lambda z w')
poly_u = v**2 + s.sqrt(2)*v*x + x**2/2
poly_v = lam*(s.sqrt(2)*v*y**3 + y**4/4)
covariance = {(x,x): xx, (y,y): yy, (x,y): xy, (y,x): xy}

def wick(labels):
    if not labels:
        return s.Integer(1)
    return sum(covariance[labels[0], labels[j]]*wick(labels[1:j]+labels[j+1:])
               for j in range(1,len(labels)))

mgf = s.exp(xx*z**2/2 + xy*z*w + yy*w**2/2)

def moment(poly, method):
    return s.expand(sum(c*method(i,j) for (i,j),c in s.Poly(poly,x,y).terms()))

pairing = lambda i,j: wick((x,)*i+(y,)*j)
derivative = lambda i,j: s.diff(mgf,z,i,w,j).subs({z:0,w:0})
for poly in (poly_u, poly_v, poly_u*poly_v):
    assert s.expand(moment(poly,pairing)-moment(poly,derivative)) == 0
u, correction, uv = (moment(p,pairing) for p in (poly_u,poly_v,poly_u*poly_v))
assert s.expand(-correction+3*lam*yy**2/4) == 0
assert s.expand(-uv+u*correction+3*lam*xy*yy*(xy+4*v**2)/2) == 0
print('R048: d_s log Z0 =', -correction)
print('R048: d_s (ZU/Z0) =', s.factor(-uv+u*correction))
print('R048: d_s ZU / Z0(0) =', s.factor(-uv))

mp.mp.dps = 85
q = mp.exp(-mp.pi)
log_product = 96*mp.fsum(mp.log1p(q**(2*k+1)) for k in range(28))
log_series = 96*mp.fsum((-1)**(m+1)*q**m/(m*(1-q**(2*m))) for m in range(1,55))
product_tail = 96*q**57/(1-q**2)
series_tail = 96*q**55/(55*(1-q)*(1-q**110))
assert abs(log_product-log_series) < product_tail+series_tail
for alpha in (mp.mpf(1)/4,mp.mpf(3)/4):
    p = sorted(abs(n+alpha) for n in range(-28,28))
    assert p[:28] == [mp.mpf(2*k+1)/4 for k in range(28)]
print('S reference trace:', mp.nstr(mp.exp(log_product),37))
print('S trace ratio floor:', mp.nstr(mp.exp(-log_product),37))
print('log tails:', mp.nstr(product_tail,12), mp.nstr(series_tail,12))
