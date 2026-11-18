"""Sourced complex scalaron Gaussian on the existing Euclidean residue branch.
No independent boson measure or new physical model is selected.
"""
from pathlib import Path
import hashlib
import json
import sympy as s

v, alpha = s.symbols('v alpha', positive=True)
F, R, xi, theta = s.symbols('F R xi theta', real=True)
mu = F + 2*alpha*R
rho = F+s.I*xi
base = s.exp(-v*(F*R+alpha*R**2))
source_ratio = s.exp(theta*mu-alpha*theta**2/v)
checks=[]
def need(name, result):
    checks.append({'name':name, 'passed':bool(result)})
    if not result:
        raise AssertionError(name)

# Direct real Gaussian Fourier integration after scaling the real variable.
x, t = s.symbols('x t', real=True)
gauss = s.integrate(s.exp(-x*x), (x, -s.oo, s.oo))
need('absolute_Gaussian_factor', gauss == s.sqrt(s.pi))
gauss_fourier = s.integrate(s.exp(-x*x+s.I*t*x), (x, -s.oo, s.oo))
need('independent_real_Gaussian_Fourier_integral', s.simplify(gauss_fourier-s.sqrt(s.pi)*s.exp(-t*t/4)) == 0)
direct = s.exp(-v*F*R+theta*F)*s.exp(-alpha*(theta-v*R)**2/v)
need('whole_sourced_identity_with_absolute_normalization', s.simplify(direct/(base*source_ratio)) == 1)
need('source_square_completion', s.expand(-v*xi**2/(4*alpha)+s.I*(theta-v*R)*xi +v*(xi-2*s.I*alpha*(theta-v*R)/v)**2/(4*alpha)+alpha*(theta-v*R)**2/v)==0)

# Independent original-R differentiation gives each direct inserted integral.
for n in range(5):
    via_xi = sum(s.binomial(n,j)*F**(n-j)*(-1/v)**j*s.diff(s.exp(-v*alpha*R**2),R,j) for j in range(n+1))/s.exp(-v*alpha*R**2)
    via_source = s.diff(source_ratio,theta,n).subs(theta,0)
    need(f'direct_auxiliary_moment_{n}_matches_sourced_integral',s.simplify(via_xi-via_source)==0)

need('first_moment_original_composite',s.diff(source_ratio,theta).subs(theta,0)==mu)
second=s.simplify(s.diff(source_ratio,theta,2).subs(theta,0))
need('neutral_second_moment_retains_negative_contact',s.simplify(second-(mu**2-2*alpha/v))==0)
need('neutral_contact_corrected_word_is_original_composite_square',s.simplify(second+2*alpha/v-mu**2)==0)
need('connected_contact',s.diff(s.log(source_ratio),theta,2)==-2*alpha/v)
need('auxiliary_sourced_Schwinger_Dyson',s.simplify(s.diff(source_ratio,theta)-(mu-2*alpha*theta/v)*source_ratio)==0)

# Raw F(H) in the actual printed homogeneous residue branch.
chi0,u=s.symbols('chi0 u',real=True)
F_orig=chi0**2/(4*s.pi**2)+u/(96*s.pi**2)
need('original_F_substitution_preserved', s.expand(mu.subs(F,F_orig)-(chi0**2/(4*s.pi**2)+u/(96*s.pi**2)+2*alpha*R))==0)

root=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
sources=[
 ('experiments/theory-contracts/craig-original-source-bridge-20261103/RESIDUE-REFERENCE-AND-FULL-INTEGRATION.txt',[[20,60],[143,173]]),
 ('experiments/theory-contracts/original-local-source-variable-background-20261102/REFERENCE-SOURCE-COMPLETION.txt',[[82,140]]),
 ('experiments/theory-contracts/original-local-source-variable-background-20261102/FULL-FUNCTOR-REVIEW.txt',[[65,125],[130,162]]),
 ('experiments/theory-contracts/original-local-source-variable-background-20261102/WEIGHT-FULL.txt',[[45,95]]),
 ('documents/archive--tfpt-45--source_extracts--04_qft_source/sections/007-subsection-renormalized-physical-and-scheme-readout-laye.tex',[[1,28]]),
]
out={'checks':checks,'passed':sum(x['passed'] for x in checks),'total':len(checks),
 'Atlas_revision_at_start':'aa7a2fae187f8d99d93468a1',
 'sourced_ratio':'exp(theta*(F+2*alpha*R)-alpha*theta**2/v)',
 'first':'F+2*alpha*R','second':'(F+2*alpha*R)**2-2*alpha/v',
 'neutral_contact':'-2*alpha/v',
 'integrated_variable':'introduced exact complex auxiliary leg only',
 'original_H_or_geometry_integral_evaluated':False,
 'new_measure_selected':False,'canonical_edits':False,
 'sources':[{'path':p,'sha256':hashlib.sha256((root/p).read_bytes()).hexdigest(),'lines':lines} for p,lines in sources]}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
