"""Event-index and physical-clock signatures of the opposite-gallery jump law."""
from pathlib import Path
import sys
import sympy as s
from flaggen import transfers
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"fortsetzung"))
from apartment import FLAGS, linear_maps

def main():
    mats=transfers()
    T=sum(mats.values(),s.zeros(21))/8
    one=s.ones(21,1)
    assert T==T.T and T*one==one
    assert T.eigenvals()=={1:1,-s.S.One/8:8,s.sqrt(2)/4:6,-s.sqrt(2)/4:6}
    v=(mats[1]-mats[2])*one/8
    corr=lambda k:(v.T*T**k*v)[0]/21
    assert corr(0)==s.Rational(1,28) and corr(1)==-s.Rational(1,168)
    assert T*T*v==v/8
    elapsed=s.symbols('t', real=True)
    a=s.sqrt(2)/4
    form=s.exp(-elapsed)*(s.cosh(a*elapsed)/28 - s.sinh(a*elapsed)/(168*a))
    assert form.subs(elapsed,0)==corr(0)
    assert s.simplify(s.diff(form,elapsed).subs(elapsed,0)/corr(0))==-s.Rational(7,6)
    assert s.simplify(s.Rational(1,28)-1/(168*a))>0
    # GL3(F2) is transitive on all ordered allowed transitions.
    seed=next((i,j) for i in range(21) for j in range(21) if T[i,j])
    orbit=set()
    for g in linear_maps():
        def rename(flag):
            p,line=flag
            return FLAGS.index((g[p],tuple(sorted(g[x] for x in line))))
        orbit.add(tuple(rename(FLAGS[i]) for i in seed))
    assert len(orbit)==168 and orbit=={(i,j) for i in range(21) for j in range(21) if T[i,j]}
    print("Ordered allowed transitions: one orbit of",len(orbit))
    print("Exact event correlations, zero to six intervening jumps:",[corr(i) for i in range(7)])
    print("Continuous unit-rate positive-time correlation:",form)
    print("Slope relative to zero-lag density:",-s.Rational(7,6))
    print("Spectral gap of the flag generator, unsigned lab time:",1-a)

if __name__=="__main__":
    main()
