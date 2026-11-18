"""Gauge-covariant local controls on the same ordered E8 source channels."""
from contextlib import redirect_stdout
from itertools import permutations, product
from pathlib import Path
import io
import runpy
import sys

root=Path(sys.argv[1])
prior=root/"experiments/theory-contracts/recursive-mirrored-source-20261103/check_seam_ordered_yukawa_20261111.py"
with redirect_stdout(io.StringIO()):
    a=runpy.run_path(str(prior))
sp=a["sp"]
corr=a["current_response"]
vertex=a["vertex_response"]
adjoint_phase=a["epsilon"]

# Color-singlet vacuum-contracted quadratic pair in the actual E8 realization.
color_roots=tuple(tail+(0,)*5 for tail in sorted(set(permutations((1,-1,0)))))
cartan=((1,-1,0)+(0,)*5,(1,1,-2)+(0,)*5)
def paired_color(roots,positions):
    out=sp.S.Zero
    for r in color_roots:
        dual=tuple(-x for x in r)
        out += adjoint_phase(dual,r)*vertex((r,dual)+roots,(4,3)+positions)
    for d in cartan:
        out += corr(roots,positions,cartan=((d,4),(d,3)))/a["charge"](d,d)
    return sp.factor(out)

assert paired_color((),())==8
x=sp.symbols("x")
cross=3*(4-x)/(4*(3-x))
rho_q=(4+2*(cross+1/cross)+sp.Rational(2,3)*(1/(4-x)-sp.Rational(1,4))*
       (1/(3-x)-sp.Rational(1,3)))/8
assert sp.series(rho_q,x,0,3)==1+x**2/432+sp.Order(x**3,x)
for x0 in (sp.Rational(1,10),sp.Rational(1,100)):
    positions=(2,sp.Rational(1,2),x0,0,-1)
    for tag,slots in a["slots"].items():
        first=a["alpha"] if tag!="u" else a["complement"][0]
        for left,right,higgs in slots:
            for i,j,k in product(range(3),repeat=3):
                p,_=a["matter_vertex"](left,i)
                q,_=a["matter_vertex"](right,j)
                h,_=(a["conjugated_down_higgs"](k) if tag=="u" else a["matter_vertex"](higgs,k))
                roots=(first,a["beta"],p,q,h)
                raw=vertex(roots,positions)
                if raw==0:
                    continue
                expected=1 if tag=="e" else rho_q.subs(x,x0)
                assert paired_color(roots,positions)/(8*raw)==expected
                assert paired_color((first,a["beta"],q,p,h),positions)/(
                    8*vertex((first,a["beta"],q,p,h),positions))==expected
                break
            else:
                raise RuntimeError(("no charged word",tag,left,right))
print("Color singlet pair: factor 1 for electron, and 1+delta^2/432+O(delta^3) for all quark colors.")
print("The same factor acts on both exchanged words. No regular local Yukawa correction survives.")

# A compact smooth Ward phase is removable in the local neutral operator
# if it is the endpoint phase of a flat electromagnetic connection.
phase,charge=sp.symbols("f q",real=True)
delta=sp.symbols("delta",positive=True)
forward=sp.exp(sp.I*phase*charge*delta)
backward=sp.exp(-sp.I*phase*charge*delta)
transport=sp.exp(-sp.I*phase*charge*delta)
assert sp.simplify(forward*transport)==1
assert sp.simplify(backward/transport)==1
pole,regular=sp.symbols("P S")
unlinked=((pole/delta+regular)*forward+(-pole/delta+regular)*backward)/2
linked=((pole/delta+regular)*forward*transport+(-pole/delta+regular)*backward/transport)/2
assert sp.limit(unlinked,delta,0)==regular+sp.I*phase*charge*pole
assert sp.simplify(linked)==regular
print("Flat compact Ward phase produces i*q*f*P only before endpoint parallel transport; the linked local coefficient is S.")

# A compact phase removed by common transport leaves E_11 decoupled at
# the quarter cross ratio. Compare it with the third diagonal entry.
mp=a["mp"]
lam=a["lY"]
phi=a["phi0"]
lepton_extreme=mp.mpf(96)/49*phi**3
largest_possible=lam**5/4
needed=lepton_extreme/largest_possible
assert needed>2
print("Every pure diagonal length split after flat phase cancellation obeys s_e/s_tau <= lambda^5/4.")
print("A family-diagonal kinetic correction needs z_tau/z_e >=",mp.nstr(needed,14))

# Can the same family currents supply that norm on native charged states?
# Zero-mode action on grade-one spinor roots is the simple-pole coefficient.
alpha,beta=a["alpha"],a["beta"]
gamma=a["complement"][0]
weight_a,weight_g=sp.symbols("a g",real=True)
def grade_one_generator(slot,terms):
    p0,_=a["matter_vertex"](slot,0)
    basis=[a["matter_vertex"](slot,j)[0] for j in range(3)]
    basis.append(p0[:5]+(sp.Rational(-1,2),)*3)
    result=sp.zeros(4)
    for j,p in enumerate(basis):
        for r,c in terms:
            if a["charge"](r,p)==-1:
                target=tuple(x+y for x,y in zip(r,p))
                result[basis.index(target),j]+=c*adjoint_phase(r,p)
    return result

terms_a=tuple((r,weight_a*c) for r,c in a["hermitian_roots"](alpha))+tuple(
    (r,weight_g*c) for r,c in a["hermitian_roots"](gamma))
terms_b=a["hermitian_roots"](beta)
all_matter=sorted({i for rows in a["slots"].values() for left,right,h in rows for i in (left,right)})
for slot in all_matter:
    A=grade_one_generator(slot,terms_a)
    B=grade_one_generator(slot,terms_b)
    assert A.T==A and B.T==B
    single=(A.T*A)[:3,:3]
    assert single==sp.diag(weight_a**2,weight_a**2,weight_g**2)
    late=A*B
    early=B*A
    assert (late.T*late)[:3,:3]==sp.diag(0,weight_g**2,weight_a**2)
    assert (early.T*early)[:3,:3]==sp.diag(weight_a**2,0,0)
    assert ((late+early).T*(late+early))[:3,:3]==sp.diag(weight_a**2,weight_g**2,weight_a**2)
    assert ((late-early).T*(late-early))[:3,:3]==sp.diag(weight_a**2,weight_g**2,weight_a**2)
print("Grade-one reflected Gram for A alone: (a²,a²,g²) on every actual matter leg.")
print("Both equal-clock orderings AB and BA give combined Gram (a²,g²,a²), even coherently or by commutator.")
print("This grade-one equal-time norm has z_tau/z_e=1 and misses the exact lepton lower bound.")

# The same invariant color pair on the external two-point Gram survives
# at a finite separation and supplies a sector scalar without a new field.
assert sp.cancel(rho_q-(1+cross+1/cross)/3)==0
TAB=a["TAB"]
rounding=a["half_printing"]
ed_lo=mp.fprod(mp.mpf(TAB[f][1])-rounding for f in ("e","mu","tau"))/mp.fprod(
    mp.mpf(TAB[f][1])+rounding for f in ("d","s","b"))
ed_hi=mp.fprod(mp.mpf(TAB[f][1])+rounding for f in ("e","mu","tau"))/mp.fprod(
    mp.mpf(TAB[f][1])-rounding for f in ("d","s","b"))
gram_required=(mp.root(ed_lo,3),mp.root(ed_hi,3))
assert gram_required[0]>1
mark_crosses=set()
for first,second,left,right in permutations((sp.S.One,sp.I,-sp.S.One,-sp.I)):
    mark_crosses.add(sp.simplify((first-left)*(second-right)/((first-right)*(second-left))))
assert mark_crosses=={sp.Rational(1,2),sp.Integer(2),sp.Integer(-1)}
assert mp.mpf(7)/6>gram_required[1]
print("External color pair gives exact Gram scalar (1+F+1/F)/3 at finite separation.")
print("Target E/D requires identical endpoint factors in",tuple(mp.nstr(r,12) for r in gram_required))
print("The positive quarter-mark cross ratios 2 and 1/2 give 7/6, outside that interval.")

# Old Ramond propagation alone, even if transferred directly to the equal
# endpoint grade-one norm model, supplies too small a relative weight.
old_d=mp.mpf(1)/20
time_weights=[mp.exp(2*(32-(30+2*quarter))*old_d) for quarter in (mp.mpf(1)/4,mp.mpf(3)/4)]
assert max(time_weights)<needed
print("The old Ramond window propagation gives weight ratios",tuple(mp.nstr(x,12) for x in time_weights),
      "as a separate conditional norm test.")
