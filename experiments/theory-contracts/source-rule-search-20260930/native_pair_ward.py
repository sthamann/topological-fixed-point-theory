"""Independently check the pair norm using original affine Ward brackets."""
from collections import defaultdict
from fractions import Fraction
import json
from pathlib import Path
import sys
import sympy as s
sys.path.insert(0,'/Users/stefanhamann/Projekte/tfpt-theoryv4')
from tfpt_explorer.current_block_geometry import _actual_current_source, _kappa, _bracket
from tfpt_explorer.neutral_source_response import _native_neutral_branch

source=_actual_current_source()
alg=source['chevalley']
roots=sorted(_native_neutral_branch()['retained'])
pairs=[]
for r in roots:
    a={alg.ridx[r]:Fraction(1)}
    opp={alg.opp[alg.ridx[r]]:Fraction(1)}
    k=_kappa(source,a,opp)
    assert k in (-1,1)
    adj={i:c/k for i,c in opp.items()}
    assert _kappa(source,a,adj)==1
    pairs.append((a,adj))

def terms(currents, labels):
    if not currents:
        return {():Fraction(1)}
    if len(currents)==1:
        return {}
    result=defaultdict(Fraction)
    for j in range(1,len(currents)):
        rest=currents[1:j]+currents[j+1:]
        rest_labels=labels[1:j]+labels[j+1:]
        k=_kappa(source,currents[0],currents[j])
        if k:
            for signature,c in terms(rest,rest_labels).items():
                result[signature+((labels[0],labels[j],2),)]+=k*c
        bracket=_bracket(source,currents[0],currents[j])
        if bracket:
            replaced=list(currents[1:]);replaced[j-1]=bracket
            for signature,c in terms(tuple(replaced),labels[1:]).items():
                result[signature+((labels[0],labels[j],1),)]+=c
    return dict(result)

total=defaultdict(Fraction)
for a,ad in pairs:
    for b,bd in pairs:
        for left in ((a,ad),(ad,a)):
            for right in ((b,bd),(bd,b)):
                for sig,c in terms(left+right,(0,1,2,3)).items():
                    total[sig]+=c/4

t,u=s.symbols('t u')
z=(t,s.Integer(0),1+u,s.Integer(1))
F=s.factor(sum((s.Rational(c.numerator,c.denominator)/s.prod((z[i]-z[j])**power for i,j,power in sig)
                for sig,c in total.items() if c),s.Integer(0)))
# Pair connected subtraction removes the pair identity contribution only.
connected=s.factor(F-s.Rational(48**2)/(t*t*u*u))
regular=s.factor(s.limit(s.limit(connected,t,0),u,0))
assert regular==624
result={'method':'original Chevalley kappa/bracket affine Ward recursion',
 'roots':48,'adjoint_even_four_current_words':len(pairs)**2*4,
 'connected_pair_function':str(connected),
 'regular_pair_pair_coefficient':str(regular),
 'primary_pair_pair_coefficient_after_subtracting_stress':'48',
 'verdict':'unprojected native connected response is 624, not 48; no physical port chosen'}
Path(__file__).with_name('native_pair_ward_results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(result,ensure_ascii=False,indent=2))
