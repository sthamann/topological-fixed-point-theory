"""Distinguish the old scalar-square theta lift from the new F20 lift."""
from pathlib import Path
import hashlib
import json
import sympy as s

repo = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
old = repo/'experiments/theory-contracts/elliptic-moebius-source-bridge-20261026/MARKED-Q-DIVISOR-CENTRAL-TWIST.txt'
pin = 'da1f75831a46a0c0662f931486b77e6e2085126df452b009d1a90351d413d904'
checks = []


def need(value, label):
    if not bool(value):
        raise RuntimeError(label)
    checks.append(label)


need(hashlib.sha256(old.read_bytes()).hexdigest() == pin, 'old marked theta source unchanged')
I = s.I
S = s.diag(1,I)
Y = s.Matrix([[0,-I],[I,0]])
Z = s.diag(1,-1)
cx = s.Matrix([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
G = s.kronecker_product(Y,Z)*cx*s.kronecker_product(s.eye(2),S)*s.kronecker_product(S,s.eye(2))
z8 = (1+I)/s.sqrt(2)
Ghat = s.simplify(G/z8)
need(G**4 == s.eye(4), 'native F20 G fourth power')
need(G**2 != s.eye(4) and G**2 != -s.eye(4), 'native F20 G has nonscalar square')
need(s.trace(G**2) == 0, 'native F20 square trace zero')
need(s.simplify(Ghat**2 + I*G**2) == s.zeros(4), 'symplectic lift square minus i times G square')
need(s.simplify(Ghat**4 + s.eye(4)) == s.zeros(4), 'both lifts can share fourth power minus I')
need(s.simplify(Ghat**2 - I*s.eye(4)) != s.zeros(4), 'symplectic F20 lift not the old theta scalar-square lift')
result = {
    'checks':checks,
    'checks_pass':True,
    'old_theta_relation_source':str(old.relative_to(repo)),
    'old_theta_relation_sha256':pin,
    'old_theta_relation':'alpha^2=+i I4; projective order2',
    'new_f20_relation':'Ghat^2=-i G^2 nonscalar; projective order4',
    'trace_of_new_square':'0',
    'trace_of_old_square':'4i',
    'verdict':'NO_CONJUGACY_EVEN_UP_TO_CENTRAL_PHASE_BETWEEN_THE_TWO_LIFTS',
    'scope':'The scalar mu8 resource is shared. These distinct group operations are not identified by equal order8 and fourth power minusI. This does not contradict either separate source construction or the native F20 representation bridge.'
}
text=json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+'\n'
Path(__file__).with_name('distinct-quarter-lifts-results.json').write_text(text)
print(text,end='')
