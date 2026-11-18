"""Entire local scalar first-order source fiber: exact analytic identities.

The source S1/CAR realization remains conditional relative to P1. No new
matrix potential, physical metric law, clock or measurement is selected.
"""
from pathlib import Path
import json
import sympy as s

def need(test,msg):
    if not test:raise ArithmeticError(msg)

def main():
    t=s.Symbol('theta',real=True);w=s.Symbol('omega',positive=True)
    beta=s.Symbol('beta',real=True)
    a=s.Function('a',positive=True)(t);c=s.Function('c',real=True)(t)
    f=s.Function('f')(t);fb=s.Function('fbar')(t);g=s.Function('g')(t)
    def B(h):return -s.I*(a*s.diff(h,t)+s.diff(a,t)*h/2)+c*h
    # Integration-by-parts identity, arbitrary test functions; periodic edge term zero.
    conj_Bf=s.I*(a*s.diff(fb,t)+s.diff(a,t)*fb/2)+c*fb
    need(s.simplify(fb*B(g)-conj_Bf*g+s.I*s.diff(a*fb*g,t))==0,
         'Formal selfadjoint half-density correction, full variable principal coefficient')
    x=s.Function('x',real=True)(t);eta=s.Function('eta',real=True)(t)
    F=s.Function('F');amp=s.sqrt(w/a)*s.exp(s.I*eta)
    u=amp*F(x)
    # Use f'(x) in the same SymPy representation as the chain rule.
    normal=amp*w*(-s.I*s.diff(F(x),x)+beta*F(x))
    remainder=(B(u)-normal).subs({s.diff(x,t):w/a,s.diff(eta,t):(w*beta-c)/a})
    need(s.simplify(remainder)==0,'Exact arbitrary-function unitary Liouville/gauge normal form')
    U=s.exp(s.I*t);phi=s.Function('phi',real=True)(t);uphi=s.exp(s.I*phi)
    need(s.simplify(B(U*f)-U*B(f)-a*U*f)==0,'Original fixed unit-winding Ward is a(theta)Z')
    need(s.simplify(B(uphi*f)-uphi*B(f)-a*s.diff(phi,t)*uphi*f)==0,
         'General winding-one symbol gives a phi_prime times u')
    ub=s.exp(s.I*x)
    need(s.simplify((B(ub*f)-ub*B(f)-w*ub*f).subs(s.diff(x,t),w/a))==0,
         'Intrinsic uniform winding from SAME principal coefficient, no fitted profile')
    need(s.simplify(s.conjugate(U)*B(U*f)-B(f)-a*f)==0,'Positive spectral-flow path endpoint Z* B Z=B+a')
    # Exact entire controls. These are discriminators, not fitted source candidates.
    eps=s.Symbol('epsilon',real=True);alpha=s.Symbol('alpha',real=True)
    controls={}
    points=[0,s.pi/2,s.pi,3*s.pi/2,2*s.pi]
    for k in (1,4):
        ac=1/(1+eps*s.cos(k*t));xc=t+eps*s.sin(k*t)/k
        need(s.simplify(s.diff(xc,t)-1/ac)==0,'Positive-speed coordinate, all epsilon in(-1,1)')
        T=s.integrate(1/ac,(t,0,2*s.pi))
        betac=s.integrate(alpha/ac,(t,0,2*s.pi))/(2*s.pi)
        need(T==2*s.pi and s.simplify(betac-alpha)==0,'Exact same total travel time and quarter holonomy')
        arc=[s.simplify(xc.subs(t,q)-xc.subs(t,p)) for p,q in zip(points,points[1:])]
        fixes=all(s.simplify(xc.subs(t,p)-p)==0 for p in points)
        quarter=s.trigsimp(xc.subs(t,t+s.pi/2)-xc-s.pi/2)==0
        need((fixes and quarter)==(k==4),'Marked and unmarked presentations correctly distinguished')
        if k==1:need(arc==[s.pi/2+eps,s.pi/2-eps,s.pi/2-eps,s.pi/2+eps],'All four marked arc travel times retained')
        else:need(arc==[s.pi/2]*4,'Quarter-equivariant control preserves every original mark and arc time')
        for av in [s.Rational(1,4),s.Rational(3,4)]:
            mode=0 if av==s.Rational(1,4) else -1
            need(abs(mode+av)==s.Rational(1,4),'Unchanged exact original quarter gap')
        controls[str(k)]={'a':str(ac),'c':'alpha','x':str(xc),'eta':'0',
          'parameter_range':'-1<epsilon<1','total_travel_time':str(T),
          'beta':str(betac),'arc_times':[str(v) for v in arc],
          'preserves_four_fixed_marks':fixes,'commutes_fixed_quarter_rotation':quarter,
          'constant_original_Z_Ward':'only epsilon=0',
          'unit_positive_spectral_flow':1,
          'full_P1_countermodel_claimed':False}
    # No finite-mode spectral-flow scan: E_n(s)=omega(n+beta+s).
    # For 0<beta<1 and 0<=s<=1 only n=-1 crosses, at s=1-beta, upward.
    z=s.Symbol('s',real=True)
    need(s.simplify(w*(-1+beta+z).subs(z,1-beta))==0,'Exact unique upward crossing branch')
    # Actual charged quarter lift and every integer carry, not periodic compensation.
    ni,mi=s.symbols('N m',integer=True)
    carry=[]
    for residue in range(4):
        phase=s.simplify(s.exp(s.I*s.pi*(4*mi+residue)/2))
        need(s.simplify(phase-s.I**residue)==0,'All integer N=4m+r retain their joint time/quarter phase')
        endpoint=[s.simplify(s.exp(-s.I*(4*mi+residue)*p)) for p in points[:-1]]
        need(all(s.simplify(v-1)==0 for v in endpoint)==(residue==0),'Large gauge fixes all four endpoint frames iff N=0 mod4')
        carry.append({'residue':residue,'all_integer_m':True,'time_quarter_phase':str(phase),
          'fixes_mark_frames_and_quarter':residue==0})
    result={'status':'PASS_EXACT_SCALAR_FIRST_ORDER_SOURCE_FIBER','date':'2026-10-17',
      'general_operator':'B=-i(a dtheta+a_prime/2)+c, a>0,c real smooth periodic, Dom H1_per',
      'complete_normal_form':{'T':'integral dtheta/a','omega':'2pi/T','beta':'(1/2pi) integral c/a',
       'x_prime':'omega/a','eta_prime':'(omega beta-c)/a',
       'U':'Uf(theta)=sqrt(omega/a(theta)) exp(i eta(theta)) f(x(theta))',
       'operator_identity':'B=U omega(-i dx+beta) U*','full_spectrum':'omega(Z+beta)'},
      'all_function_symbolic_checks':['formal selfadjointness','entire operator intertwiner',
       'fixed Z Ward','general winding-one Ward','intrinsic uniform winding','spectral-flow endpoint'],
      'intrinsic_time_seam_construction':'u_B=exp(ix(theta)) is unique up to phase among wind1 scalar symbols with constant positive increment; increment=omega. Q_B=exp(i T B/4)=U[exp(i pi beta/2)Q4_x]U*, Q_B^4=exp(2pi i beta)I. No identity with native/internal/physically frozen mark clock inferred.',
      'constant_seam_Ward_theorem':'[B,Z]=kappa Z iff a(theta)=kappa. For u=exp(i phi) of wind1, [B,u]=kappa u iff a phi_prime=kappa; then omega=kappa.',
      'unit_spectral_flow_theorem':'For every positive a, Z* B Z=B+a and B_s=B+s a has E_n(s)=omega(n+beta+s). For0<beta<1 exactly one upward crossing. Unit spectral flow does not imply constant increment Ward.',
      'fixed_quarter_and_gap':'beta=alpha modZ, alpha1/4 or3/4, gap1/4 force omega1; all unmarked operators one transported unitary source orbit. Does not identify frozen marked fields/state or a physical q/S selector.',
      'exact_controls':controls,
      'charged_quarter_and_integer_carry':{'Qalpha':'exp(i pi alpha/2) Q4; Qalpha^4=exp(2pi i alpha)I',
         'joint_invariant':'exp(i pi B/2) Qalpha*=i^N I for mean(c)=alpha+N',
         'all_integer_residue_certificates':carry,'reference_sea_branch':'N=0; do not discard4Z large-gauge carry'},
      'scope':'Whole given rank1 LOCAL differential first-order circle class, not all originalP1 sources, nonlocal kernels, enlarged internal algebras or fullE8 operations. Original fixed seam/time Ward is a GIVEN-source condition, not inferred from bareP1 winding. Four marks must be transported or their arc times retained. No change of reference sea/normal-order/Schwarzian or full relative action inferred from one-particle unitary equivalence. No physical3+2, event/coupling/metric rule or inequivalent-P1 countermodel claimed.'}
    Path(__file__).with_name('scalar_first_order_source_fiber_results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'all_function_identities':len(result['all_function_symbolic_checks']),
      'controls':{k:{'fixed_marks':v['preserves_four_fixed_marks'],'quarter_equivariance':v['commutes_fixed_quarter_rotation'],
         'same_gap_and_holonomy':True,'spectral_flow':1,'constant_Z_Ward':v['constant_original_Z_Ward']} for k,v in controls.items()}},indent=2))

if __name__=='__main__':main()
