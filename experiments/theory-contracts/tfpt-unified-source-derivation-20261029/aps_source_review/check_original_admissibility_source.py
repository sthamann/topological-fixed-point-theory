"""Actual four-generator admissibility complex, not a chosen spacetime model."""
from pathlib import Path
import json
import sympy as s

checks=[]
def need(name, value):
    good=bool(value)
    checks.append({'name':name,'passed':good})
    if not good: raise AssertionError(name)

n=4; dim=1<<n; I=s.eye(dim); zero=s.zeros(dim)
cre=[]
for j in range(n):
    a=s.zeros(dim)
    for mask in range(dim):
        if not (mask>>j)&1:
            a[mask|(1<<j),mask]=(-1)**((mask&((1<<j)-1)).bit_count())
    cre.append(a)
ann=[a.T for a in cre]
need('all creation anticommutators',all(cre[i]*cre[j]+cre[j]*cre[i]==zero for i in range(n) for j in range(n)))
need('all mixed CAR anticommutators',all(ann[i]*cre[j]+cre[j]*ann[i]==(I if i==j else zero) for i in range(n) for j in range(n)))
k=s.symbols('k0:4',real=True); z,u=s.symbols('z u')
Q=sum((k[j]*cre[j] for j in range(n)),zero); Qt=Q.T
K=sum(a*a for a in k); D=Q+Qt
need('whole four-penalty differential is nilpotent',Q*Q==zero)
need('whole Laplacian equals sum squares on every ghost grade',Q*Qt+Qt*Q==K*I)
need('actual selfadjoint odd source square',D*D==K*I)
need('contracting homotopy numerator on every positive joint fiber',Q*Qt+Qt*Q==K*I)
need('zero-penalty fiber keeps all sixteen ghosts',D.subs(dict.fromkeys(k,0))==zero)
grades=[mask.bit_count() for mask in range(dim)]
mult=[grades.count(j) for j in range(n+1)]
need('actual grade multiplicities 1 4 6 4 1',mult==[1,4,6,4,1])
need('whole degree character',(s.expand(sum(u**a for a in grades)-(1+u)**4))==0)
need('index supertrace vanishes',sum((-1)**a for a in grades)==0)
need('degree weighted supertrace vanishes',sum((-1)**a*a for a in grades)==0)
need('two equal parity sectors',sum(a%2==0 for a in grades)==8 and sum(a%2==1 for a in grades)==8)
need('whole energy resolvent numerator',(D-z*I)*(D+z*I)==(K-z*z)*I)
need('whole degree-weighted resolvent trace',s.expand(sum(u**grades[a]*(D+z*I)[a,a] for a in range(dim))-z*(1+u)**4)==0)
need('all original constraint insertions survive in off-diagonal source',all((D+z*I)[1<<j,0]==k[j] for j in range(n)))
need('ghost grade zero resolvent diagonal retains energy denominator',(D+z*I)[0,0]==z)
need('ordinary source trace has no unsquared constraint direction',s.trace(D)==0)

# Exact joint spectral projection calculation: every commuting projection pattern.
lam=s.symbols('lambda0:4',positive=True)
patterns=[tuple((mask>>j)&1 for j in range(n)) for mask in range(dim)]
costs=[sum(lam[j]*e[j] for j in range(n)) for e in patterns]
need('only joint zero pattern is harmonic',all((cost==0)==(mask==0) for mask,cost in enumerate(costs)))
need('every positive pattern is bounded below by a realized positive penalty',all(cost==0 or any(s.simplify(cost-lam[j]).is_nonnegative for j in range(n) if patterns[mask][j]) for mask,cost in enumerate(costs)))
need('at most sixteen fixed finite penalty eigenvalues',len(costs)==16 and all(s.simplify(sum(lam)-cost).is_nonnegative for cost in costs))
need('positive penalty reweighting preserves the whole joint kernel',all((sum((j+2)*lam[j]*e[j] for j in range(n))==0)==(mask==0) for mask,e in enumerate(patterns)))

# Resolvent source variation identity, with independent variations in every K_j.
v=s.symbols('v0:4'); dD=sum((v[j]*(cre[j]+ann[j]) for j in range(n)),zero)
dK=2*sum(k[j]*v[j] for j in range(n))
need('full energy-source variation including off-diagonal insertions',s.simplify((dD*(K-z*z)-(D+z*I)*dK)+(D+z*I)*dD*(D+z*I))==zero)

out={'date':'2026-10-31','results':['R049','R050'],'scope':'Actual original four-ghost Koszul admissibility operator; no physical ghost measure or normal-sphere to spacetime identification. Symbolic scalar joint fibers encode the whole commuting spectral family.','all_passed':all(x['passed'] for x in checks),'checks':checks,'grade_multiplicities':mult}
Path(__file__).with_name('ORIGINAL-ADMISSIBILITY-SOURCE.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'all_passed':out['all_passed'],'checks':len(checks),'grade_multiplicities':mult}))
