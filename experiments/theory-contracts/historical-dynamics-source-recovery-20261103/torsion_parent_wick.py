"""Continue the supplied P25 auxiliary action, keeping actual Clifford signs.

Local flat chart only. A finite regulator proves the Gaussian identity;
it does not supply the whole current TFPT measure or curved Wick domain.
"""
from pathlib import Path
import ast, hashlib, json
import sympy as s

owner=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4/experiments/theory-contracts/craig-original-source-bridge-20261103/actual_product_pairing.py')
nodes=[]
for n in ast.parse(owner.read_text()).body:
    if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and isinstance(n.value.func,ast.Name) and n.value.func.id=='ck':
        continue
    nodes.append(n)
    if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='g5' for t in n.targets):
        break
ns={'__file__':str(owner)}
exec(compile(ast.Module(body=nodes,type_ignores=[]),str(owner),'exec'),ns)
g,g5=ns['gam'],ns['g5']
checks={}
def ck(name,x):
    ok=all(s.simplify(v)==0 for v in x) if isinstance(x,s.MatrixBase) else s.simplify(x)==0
    checks[name]=bool(ok)
    if not ok: raise AssertionError(name)

gL=[g[3],s.I*g[0],s.I*g[1],s.I*g[2]]
g5L=s.I*gL[0]*gL[1]*gL[2]*gL[3]
ck('Lorentz_chirality_is_negative_tangent_first_Euclidean',g5L+g5)
SL=s.symbols('S0:4',real=True)
SE=[-s.I*SL[1],-s.I*SL[2],-s.I*SL[3],SL[0]]
lor=s.zeros(4);eu=s.zeros(4)
emap=[1,2,3,0] # E index 1,2,3,4 -> L index 1,2,3,0
for mu in range(4):
    for a in range(4):
        for b in range(4):
            kl=sum(s.LeviCivita(mu,a,b,v)*SL[v]/2 for v in range(4))
            lor+=s.I*gL[mu]*kl*(gL[a]*gL[b]-gL[b]*gL[a])/8
            # g_L^continued=-g_E. Lowering the connection difference
            # supplies an additional minus beyond the 3-form pullback.
            ml,al,bl=emap[mu],emap[a],emap[b]
            pull=(-s.I)**sum(j==3 for j in (mu,a,b))
            ke=-pull*sum(s.LeviCivita(ml,al,bl,v)*SL[v]/2 for v in range(4))
            ke_dual=sum(s.LeviCivita(mu,a,b,v)*SE[v]/2 for v in range(4))
            ck(f'Wick_contorsion_Hodge_{mu}{a}{b}',ke-ke_dual)
            eu+=g[mu]*ke*(g[a]*g[b]-g[b]*g[a])/8
slashSL=SL[0]*gL[0]-sum((SL[i]*gL[i] for i in range(1,4)),s.zeros(4))
slashSE=sum((SE[i]*g[i] for i in range(4)),s.zeros(4))
ck('actual_Lorentz_spinlift',lor+s.Rational(3,4)*slashSL*g5L)
ck('actual_Euclidean_traditional_spinlift',eu+s.Rational(3,4)*slashSE*g5)
ck('fermion_Lagrangian_Wick_sign',eu+lor)
ck('Lorentz_axial_square_continues_to_Euclidean_square',SL[0]**2-sum(x*x for x in SL[1:])-sum(x*x for x in SE))
u=s.symbols('u1:5',real=True)
ck('derivative_axial_mixing_Wick',s.I*u[3]*SL[0]+sum(u[i]*SL[i+1] for i in range(3))-s.I*sum(u[i]*SE[i] for i in range(4)))

V=(s.eye(4)-s.I*g5)/s.sqrt(2)
ck('constant_fermion_rotation_mass',s.I*V*g5*V-s.eye(4))
for mu in range(4):
    ck(f'constant_fermion_rotation_principal_{mu}',s.I*V*(-s.I*g[mu])*V-g[mu])
    ck(f'constant_fermion_rotation_torsion_{mu}',s.I*V*(s.I*g[mu]*g5)*V+g[mu]*g5)

A,B,c0=s.symbols('A B c0',positive=True)
S,j,da=s.symbols('S J da') # J is Grassmann-even; component polynomial identity
source=s.I*B*da+c0*j
bare=A*S**2-S*source
completed=A*(S-source/(2*A))**2-source**2/(4*A)
ck('exact_auxiliary_square_with_current',bare-completed)
eff=B**2*da**2/(4*A)-s.I*B*c0*da*j/(2*A)-c0**2*j**2/(4*A)
ck('exact_effective_kinetic_current_contact',eff+source**2/(4*A))
mu,kap,lam=s.symbols('mu kappa lambda',positive=True)
subs={A:3*mu**2/(4*kap**2),B:lam*mu**2/2,c0:s.Rational(3,4)}
c=B*c0/(2*A)
F2=s.simplify((B**2/(2*A))/(4*c*c))
ck('original_physical_phase_F_squared',F2.subs(subs)-2*mu**2/(3*kap**2))
ck('original_four_fermion_contact', (c0*c0/(4*A)).subs(subs)-3*kap**2/(16*mu**2))
ck('imaginary_Euclidean_bosonic_saddle',s.simplify((s.I*B/(2*A))/(2*c)).subs(subs)-2*s.I/3)
# Integral with external source ell: (pi/A)^2 exp((ell+iB da+c0 J)^2/(4A)).
# Differentiation therefore retains every auxiliary-field insertion and contact.
ell=s.symbols('ell')
ck('original_auxiliary_insertion_contact',s.diff((source+ell)**2/(4*A),ell,2)-1/(2*A))
th=s.symbols('dtheta')
ck('effective_E_is_one_fixed_square',eff.subs(da,th/(2*c))-F2*(th-s.I*j/(2*F2))**2/2)
effL=(B*da-c0*j)**2/(4*A)
ck('effective_L_is_one_fixed_square',effL.subs(da,th/(2*c))-F2*(th-j/(2*F2))**2/2)

out={
 'result':'P25-ORIGINAL-AUXILIARY-WICK-INTEGRATION-20261104',
 'checks':checks,'passed':sum(checks.values()),
 'original_constructor_sha256':hashlib.sha256(owner.read_bytes()).hexdigest(),
 'Euclidean_orientation':'epsilon1234=+1; g5E=g1g2g3g4; g5L=-g5E',
 'metric_lowered_contorsion':'K_E=-pullback(K_L), because g_L_cont=-g_E',
 'component_map':'S_E4=S_L^0; S_Ei=-i S_L^i',
 'parent_E':'A S_E^2 - S_E.(iB da+c0 J5E)',
 'integral_per_four_vector_cell':'(pi/A)^2 exp[(ell+iB da+c0 J5E)^2/(4A)]',
 'effective_E':'B^2/(4A)(da)^2-i B c0/(2A) da.J5E-c0^2/(4A) J5E^2',
 'effective_E_single_square':'Ftheta^2/2 (dtheta-i J5E/(2 Ftheta^2))^2',
 'effective_L_single_square':'Ftheta^2/2 (dtheta-J5L/(2 Ftheta^2))^2',
 'physical_phase':'theta=2 c a; c=B c0/(2A)',
 'Ftheta_squared':'2/(3 kappa_phys^2)',
 'bosonic_saddle':'S_E*= (2i/3) dtheta, not a real geometric gradient',
 'saddle_is_full_functional_identity':False,
 'actual_full_Phi_blocks_replaced':False,
 'full_primitive_TFPT_caller_identified':False,
 'physical_Htrpos_or_spectrum_claimed':False,
 'scope':'The supplied original P25 local auxiliary polynomial and minimal spin coupling in its already supplied flat 4D chart. Wick/gamma convention derived explicitly. Exact finite-regulator Gaussian with insertions and determinant retained. Real Euclidean auxiliary contour gives the analytic Feynman Gaussian continuation of this parent; no global curved TFPT contour, complete measure, boundary or primitive source equality is inferred.'
}
Path(__file__).with_name('originale-torsion-wick-integration.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='checks'},ensure_ascii=False,indent=2))
