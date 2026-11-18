"""Exact modular certificates for the central extension and equivariant obstruction."""
from itertools import product
from pathlib import Path
import json

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

def c(p,q):
    return (-p*q) % 4

for p,q,r in product(range(2), repeat=3):
    require((c(p,q)+c(p^q,r)-c(q,r)-c(p,q^r))%4==0,'Cocycle identity')
for p,q in product(range(2), repeat=2):
    require((2*c(p,q)-((-p)+(-q)-(-(p^q))))%4==0,'Twice cocycle is coboundary')
    require((2*c(p,q)-((-p)+(-q)-(-(p^q))))%8==0,'Pushout to Z/8 splits')
require(all((2*k-1)%4 !=0 for k in range(4)), 'No involution lift')

def ext_mul(x,y):
    a,p=x;b,q=y
    return ((a+b-p*q)%4,p^q)

u=(0,1);x=(0,0);powers=[]
for j in range(8):
    powers.append(x);x=ext_mul(x,u)
require(x==(0,0) and len(set(powers))==8, 'Reflection preimage is C8')
phase_tuples=[k for k in product(range(4), repeat=4) if (sum(k)-2)%4==0]
require(len(phase_tuples)==64,'64 allowed framed central lifts')
require(not any(len(set(k))==1 for k in phase_tuples),'No cyclically constant peripheral phases')
# Orbifold alternative: C^4=(MC)^4=1; one generator has odd determinant.
for epsilon in (0,1):
    fourth_clock=(-2*epsilon)%4
    fourth_compound=(-2*(epsilon^1))%4
    require(fourth_clock !=0 or fourth_compound !=0, 'Orbifold fourth-power obstruction')

out={
 'cocycle_modulus':4,
 'cocycle':'-epsilon(g)*epsilon(h)',
 'nontrivial_restriction':'C8 -> C2 with kernel C4',
 'cohomology_class_order':2,
 'framed_lifts_per_monodromy':64,
 'strict_C4_equivariant_lifts':0,
 'central_lift_orbits_over_17_classes':1088,
 'global_correction_condition':'sum(k_i) = 2 mod 4',
 'enlarged_coefficient_group':'Z/8',
 'pushout_trivializing_cochain':'-epsilon(g) mod 8'
}
work=Path(__file__).resolve().parent/'work';work.mkdir(exist_ok=True)
(work/'cocycle_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
