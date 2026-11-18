from functools import lru_cache
import sympy as s
c=s.symbols('c',real=True)
@lru_cache(None)
def expectation(word):
    if not word: return s.Integer(1)
    a=word[0]
    if a<=0: return s.Integer(0)
    out=s.Integer(0)
    for j,b in enumerate(word[1:],1):
        before=word[1:j]; after=word[j+1:]
        out+=(a-b)*expectation(before+(a+b,)+after)
        if a+b==0:
            out+=c*s.Rational(a**3-a,12)*expectation(before+after)
    return s.expand(out)
def gram(basis):
    return s.Matrix([[expectation(tuple(reversed(lam))+tuple(-k for k in mu)) for mu in basis] for lam in basis])
# Vacuum Verma quotient by L_-1|0>=0. The remaining partitions span each grade.
assert expectation((2,-2))==c/2
assert expectation((3,-3))==2*c
assert s.expand(expectation((2,2,-2,-2))-c*(c+8)/2)==0
G=gram([(6,),(4,2),(3,3),(2,2,2)])
det=s.factor(G.det())
expected_det=s.Rational(3,4)*c**4*(2*c-1)*(5*c+22)**2*(7*c+68)
assert s.expand(det-expected_det)==0
v=s.Matrix([-108,-264,93,64])
norm=s.factor((v.T*G*v)[0])
assert s.expand(norm-12*c*(2*c-1)*(128*c+6943))==0
assert G.subs(c,s.Rational(1,2))*v==s.zeros(4,1)
assert norm.subs(c,s.Rational(1,4))<0
assert s.Rational(45,9)+s.Rational(15,5)==8
import json
from pathlib import Path
result={"verdict":"PASS_CONDITIONAL_ALGEBRA", "gram_basis":[[6],[4,2],[3,3],[2,2,2]], "gram_matrix":[[str(x) for x in row] for row in G.tolist()], "determinant":str(det), "witness":[int(x) for x in v], "witness_norm":str(norm), "excluded_interval":"0 < c_residual < 1/2: all factors except (2*c_residual-1) are positive", "sugawara_c":8, "strict_saturation_bound":"c_W < 17/2", "physical_source_computed":False, "continuum_upper_bound_proved":False}
Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
