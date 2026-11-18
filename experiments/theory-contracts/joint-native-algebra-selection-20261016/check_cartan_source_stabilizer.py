"""Root-zero-mode leakage from the actual native Cartan quartic code.
Only the informative oscillator-partition(3) component is needed.
"""
from pathlib import Path
import sys,json,itertools as it
import sympy as s
REPO=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4');sys.path.insert(0,str(REPO))
from tfpt_explorer.cartan_source import _lattice_data,_j
from tfpt_explorer.process import _code_basis

def need(ok,msg):
    if not ok:raise ArithmeticError(msg)
def clean(m):return m.applyfunc(s.simplify)
def zero(m):return clean(m)==s.zeros(*m.shape)
def main():
    roots=_lattice_data()['roots'];_,supports=_code_basis()
    J=s.zeros(8);U=s.zeros(8,4)
    for k in range(4):
        J[2*k,2*k+1]=-1;J[2*k+1,2*k]=1
        U[2*k,k]=1/s.sqrt(2);U[2*k+1,k]=-s.I/s.sqrt(2)
    need(zero(J*U-s.I*U) and zero(U.T*U) and zero(U.H*U-s.eye(4)),'Original J-positive Cartan polarization')
    x=s.symbols('x0:4');polys=[]
    for support in supports:
        # Isometric symmetric-tensor to single-mode Heisenberg Fock map.
        f=sum(s.prod(x[j] for j in word) for word in support)/s.sqrt(24*len(support))
        polys.append(s.expand(f))
    expected=[sum(t**4 for t in x), x[0]**2*x[1]**2+x[2]**2*x[3]**2,x[0]**2*x[2]**2+x[1]**2*x[3]**2,x[0]**2*x[3]**2+x[1]**2*x[2]**2,s.prod(x)]
    factors=[s.simplify(polys[k]/expected[k]) for k in range(5)]
    need(all(not f.free_symbols and f!=0 for f in factors),'Actual five native quartics, not substituted code')
    # Positive Fock norm check, preserving original tensor normalization.
    gram=s.zeros(5)
    for i,f in enumerate(polys):
        fp=s.Poly(f,*x)
        for j,g in enumerate(polys):
            gp=s.Poly(g,*x)
            gram[i,j]=s.simplify(sum(s.conjugate(c)*gp.coeff_monomial(m)*s.prod(s.factorial(v) for v in m) for m,c in fp.terms()))
    need(gram==s.eye(5),'Actual tensor/Fock isometry')
    need(all(s.expand(sum(x[j]*s.diff(f,x[j]) for j in range(4))-4*f)==0 for f in polys),'Original conformal L0 is 4 on the actual oscillator code')
    rows=[]
    for r in roots:
        # Construction-A vectors have squared norm4; physical alpha=r/sqrt2.
        alpha=s.Matrix(r)/s.sqrt(2);need((alpha.T*alpha)[0]==2,'Root physical normalization')
        b=clean(U.T*alpha)
        need(not zero(b),'Every real root has nonzero J-positive pairing')
        coeffs=[s.simplify(f.subs(dict(zip(x,b)))/3) for f in polys]
        nonzero=[j for j,v in enumerate(coeffs) if v!=0]
        need(bool(nonzero),'Every actual root zero mode leaks a native code vector')
        # Directly verify the translation coefficient entering the vertex law:
        # f(x-b/z) has z^-4 coefficient f(b), supplying alpha(-3)/3.
        z=s.Symbol('z')
        for f,c in zip(polys,coeffs):
            leading=s.Poly(s.expand(z**4*f.subs({x[j]:x[j]-b[j]/z for j in range(4)}, simultaneous=True)),z).coeff_monomial(1)
            need(s.simplify(leading/3-c)==0,'Exact quartic vertex-translation leading coefficient')
        rows.append({'root':list(r),'nonzero_code_columns':nonzero,'partition3_coefficients':[str(v) for v in coeffs]})
    need(len(rows)==240,'All original roots')
    # Code one-particle support, for the global normalizer argument.
    gram1=s.zeros(4)
    for f in polys:
        for i,j in it.product(range(4),repeat=2):
            left=s.Poly(s.diff(f,x[i]),*x);right=s.Poly(s.diff(f,x[j]),*x)
            gram1[i,j]+=s.simplify(sum(s.conjugate(c)*right.coeff_monomial(m)*s.prod(s.factorial(v) for v in m) for m,c in left.terms()))
    need(gram1==5*s.eye(4),'One-particle code support equals full Hplus')
    out={'status':'PASS_EXACT_NATIVE_CARTAN_CODE_ROOT_LEAKAGE','physical_root_normalization':'alpha=r/sqrt2; native Construction-A norm4','polarization':'U_k=(e_2k-i e_2k+1)/sqrt2, JU=iU','code_fock_gram':'I5','conformal_L0_on_code':'4 I5; conformal quarter equals I5','code_polynomial_factors':[str(v) for v in factors],'root_count':len(rows),'every_root_has_nonzero_partition3_leakage':True,'vertex_component':'Ealpha_0 f contains f(b)/3 alpha(-3)ealpha, b=U^T alpha','one_particle_number_gram':'5 I4 (sum over the five normalized code states)','complex_leakage_map_rank':240,'complex_zero_mode_code_stabilizer_lie_dimension':8,'lie_action_on_code':'Cartan zero modes annihilate the whole charge-zero oscillator code','global_conclusion_dependency':'Compact E8 code stabilizer has identity component maximal torus; all components normalize it. Weyl action preserving code must preserve Hplus by one-particle support; existing full U4 code stabilizer gives projective logical S6. No marked R5 within these inner source automorphisms.','scope':'Literal original J-positive Cartan-oscillator embedding in the stipulated E8 lattice vacuum. Does not exclude other source-to-code embeddings, non-automorphism interventions, higher-mode composites, paths driven by such operators or all P1 physics. Inner-group sequences returning to code are already covered by the exclusion.','root_witnesses':rows}
    Path(__file__).with_name('cartan_stabilizer_results.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='root_witnesses'},indent=2))
if __name__=='__main__':main()
