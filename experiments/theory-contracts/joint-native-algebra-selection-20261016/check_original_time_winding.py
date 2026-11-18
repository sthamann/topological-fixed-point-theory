"""Exact all-integer source winding/time identities. No finite clock or readout."""
import json
from pathlib import Path
import sympy as s

def need(ok, message):
    if not ok: raise RuntimeError(message)

def main():
    p=s.symbols('p0:8', real=True); alpha,q,n,m,z=s.symbols('alpha q n m z', real=True)
    Q=sum(p); N=s.Symbol('N',real=True)
    def F(x):return x*(x-1)/16+alpha*x/8
    def K(v):return N+sum(x*x/2+((k+alpha)/8-s.Rational(1,2))*x for k,x in enumerate(v))
    def c(nn,qq):return nn*qq/8+nn*(nn-1)/16+alpha*nn/8
    need(s.expand(c(m,q+n)+c(n,q)-c(n+m,q))==0,'All-integer cocycle')
    invariant=s.expand(K(p)-F(Q))
    need(s.diff(invariant,alpha)==0,'Same source invariant is independent of quarter phase')
    residues=[]
    for r in range(8):
        out=[p[(k-r)%8]+z+int(k<r) for k in range(8)]
        nn=8*z+r
        need(s.expand(sum(out)-Q-nn)==0,'Retain arbitrary integer charge carry')
        need(s.expand(K(out)-K(p)-c(nn,Q))==0,'All n=8z+r, not bounded winding scan')
        need(s.expand(K(out)-F(sum(out))-invariant)==0,'Original folded invariant')
        residues.append({'r':r,'all_integer_z':True})
    # Direct chronology versus normal form, with symbolic input charge/energy.
    e=s.Symbol('e',real=True); ts=s.symbols('t0:4',real=True); ns=s.symbols('n1:4',integer=True)
    cumulative=[s.Integer(0)]
    for v in ns:cumulative.append(cumulative[-1]+v)
    direct=sum(t*(e+c(v,q)) for t,v in zip(ts,cumulative))
    normal=sum(ts)*e+q*sum(t*v for t,v in zip(ts,cumulative))/8+sum(t*F(v) for t,v in zip(ts,cumulative))
    need(s.expand(direct-normal)==0,'Ordered time-winding normal form')
    need(s.expand(c(-n,-q).subs(alpha,1-alpha)-c(n,q))==0,'Opposite quarter and charge conjugation preserve source time')
    # Actual first two affine source charges, exact original M.
    import check_mixed_return as old
    fold=json.loads((old.R/'experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json').read_text())
    M=s.Matrix(fold['explicit_Weyl_matrix']); a1=M[:,0];a2=M[:,0]+M[:,1]
    need((a1.T*a1)[0]==1 and a2==s.ones(8,1)/2,'Actual odd vector and even root sectors')
    data={'status':'PASS_EXACT_ALL_WORD_SOURCE_TIME_WINDING',
      'date':'2026-10-17','F_alpha_q':str(F(q)),'invariant':str(invariant),
      'all_integer_winding_certificates':residues,'symbolic_chronology_checked':True,
      'opposite_quarter_transport_checked':True,'odd_first_charge':[str(x) for x in a1],
      'even_second_charge':[str(x) for x in a2],
      'normal_form':'W_N exp(-i T K) exp(-i Q sum_j(t_j N_j)/8) exp(-i sum_j t_j F_alpha(N_j))',
      'theorem_scope':'All finite words in the GIVEN global free chiral source time and W_n (and global Q phases), original CAR sea and full charge carry. Not all local fields or P1 interventions. Fixed-grade neutral cyclic source multiplicity is unchanged; raw returned field coordinates retain the supplied B_k metric. No executed compression or closed-C5 axiom. Odd W1 is not a local even E8 field; W2 restriction is on actual even CAR.'}
    Path(__file__).with_name('original_time_winding_results.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(data,indent=2))

if __name__=='__main__':main()
