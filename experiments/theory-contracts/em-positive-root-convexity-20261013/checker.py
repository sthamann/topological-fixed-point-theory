"""Symbolic identities and original interval bracket for the analytic proof."""
import hashlib
import json
from pathlib import Path
import sys
import mpmath as mp
from mpmath import iv
import sympy as s

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'verification'))
from v3_em_alpha import make_F,F_iv


def run():
    manifest=json.loads(Path(__file__).with_name('source_manifest.json').read_text())
    for source in manifest['sources']:
        if hashlib.sha256((ROOT/source['path']).read_bytes()).hexdigest()!=source['sha256']:
            raise RuntimeError('Pinned source changed: '+source['path'])
    checks=[]
    def check(name,condition):
        checks.append({'name':name,'ok':bool(condition)})
    q,p,k=s.symbols('q p k',positive=True)
    u=q*(1-q)**(-s.Rational(5,4))
    D=lambda f:s.diff(f,q)*(-2*q)
    ell=s.log(q)-s.Rational(5,4)*s.log(1-q)
    check('exact log u second derivative',s.simplify(D(D(ell))-5*q/(1-q)**2)==0)
    # Generic u(x)=exp(ell): a positive mixture with a constant.
    U,L1,L2=s.symbols('U L1 L2',positive=True)
    w=U/(p+U)
    direct=U*(L2+L1**2)/(p+U)-(U*L1/(p+U))**2
    check('log mixture convexity identity',s.simplify(direct-w*L2-w*(1-w)*L1**2)==0)
    check('native q mixture identity',s.simplify(D(D(s.log(p+u)))-
          (u/(p+u)*D(D(ell))+u*p/(p+u)**2*D(ell)**2))==0)
    check('rational d bound',s.Rational(48,24**4)==s.Rational(1,6912)<s.Rational(1,1000))
    check('rational response bound',s.Rational(1000,999)**2<2 and s.Rational(1,18)+s.Rational(1,500)<1)
    x,c,C,H=s.symbols('x c C H',positive=True)
    h=x**3-2*c**3*x**2
    check('convexity lower endpoint',s.diff(h,x,2).subs(x,2*c**3)==8*c**3)
    check('negative prethreshold polynomial factor',s.expand(x*x*(x-2*c**3))==h)
    iv.dps=60
    lo,hi=iv.mpf('0.0072973525622097'),iv.mpf('0.0072973525622100')
    flo,fhi=F_iv(lo),F_iv(hi)
    check('original F_iv strict bracket',flo.b<0 and fhi.a>0)
    mp.mp.dps=70
    native=make_F(1/(8*mp.pi),41)
    root=mp.findroot(native,mp.mpf('0.0073'))
    check('native numerical root inside certified bracket',mp.mpf('0.0072973525622097')<root<mp.mpf('0.0072973525622100'))
    check('native c3v argument unused acknowledged',native(mp.mpf('0.007'))==make_F(mp.mpf(999),41)(mp.mpf('0.007')))
    result={'verdict':'PARTIAL','mathematical_verdict':'UNIQUE_SIMPLE_POSITIVE_ROOT_PROVED_FOR_FIXED_FORMULA',
            'passed':sum(c['ok'] for c in checks),'total':len(checks),'checks':checks,
            'source_sha256':hashlib.sha256((ROOT/'verification/v3_em_alpha.py').read_bytes()).hexdigest(),
            'root_approximation':mp.nstr(root,50),'reciprocal_approximation':mp.nstr(1/root,50),
            'interval_signs':{'lower':str(flo),'upper':str(fhi)},
            'scope':'Analytic proof is in PROOF.txt; finite checks verify symbolic identities and original interval bracket. Physical inputs and root interpretation remain conditional.'}
    Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    if result['passed']!=result['total']:raise SystemExit(1)


if __name__=='__main__':run()
