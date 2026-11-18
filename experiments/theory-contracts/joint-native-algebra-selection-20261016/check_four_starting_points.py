"""Four native starts in the SAME original fold state, without a selector.

Exact source-relative response invariants, full joint quartic/erasure kernel,
actual native source transports and charged-root/time diagnostics.
"""
from pathlib import Path
import sys, json, math, itertools as it
from collections import Counter, deque
import sympy as s

R = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
BASE = R/'experiments/theory-contracts/joint-native-algebra-selection-20261016'
sys.path[:0] = [str(BASE), str(R)]
from tfpt_explorer import cartan_source as cs
import check_mixed_return as old

def need(ok, msg):
    if not ok:
        raise ArithmeticError(msg)

def cl(A):
    # Entries lie in Q(i,sqrt2,sqrt3), with only rational denominators.
    # Exact expansion reduces this ring without repeated heuristic simplify.
    return A.applyfunc(s.expand)

def mat(A):
    return [[str(v) for v in row] for row in A.tolist()]

def multi(r):
    return list(it.combinations_with_replacement(range(4), r))

def evaluate(bank, x, b):
    sub = dict(zip(x, b))
    return [cl(E.subs(sub, simultaneous=True)) for E in bank]

def polar(N, q):
    vv = [N[:,j]-(s.sqrt(6)-1)*N[:,q]/5 for j in (0,1)]
    return cl(sum((v*v.T for v in vv), s.zeros(5)))

def source_transports(N, polys, x, O):
    ld = cs._lattice_data()
    roots = ld['roots']; ri = {r:i for i,r in enumerate(roots)}
    T = s.Matrix.hstack(*map(s.Matrix,cs._TARGET_SIMPLE))
    targets = [tuple(T*s.Matrix(z)) for z in ld['root_coordinates']]
    need(set(targets)==cs._standard_doubled_roots(), 'All actual E8 charges retained')
    rays = cs._rays(roots)
    gens=[]; logs=[]
    A=[s.expand(sum(N[i,j]*polys[i] for i in range(5))) for j in range(6)]
    for ix in (0,13,2,3,1):
        ray=rays[ix]
        p=bytes(ri[cs._native_reflection(r,ray)] for r in roots)
        real=s.Matrix.hstack(*[s.Matrix(cs._native_reflection(r,ray)) for r in cs._SIMPLE])*ld['basis'].inv()
        g=s.Matrix(4,4,lambda i,j:real[2*i,2*j]+s.I*real[2*i+1,2*j])
        images=[s.expand(a.subs(dict(zip(x,g.T*s.Matrix(x))),simultaneous=True)) for a in A]
        lp=tuple(A.index(a) for a in images)
        gens.append(p);logs.append(lp)
    identity=bytes(range(240));group={identity:tuple(range(6))};queue=deque([identity])
    translations=[bytes(list(p)+list(range(240,256))) for p in gens]
    while queue:
        p=queue.popleft();lp=group[p]
        for tr,lg in zip(translations,logs):
            pp=p.translate(tr);ll=tuple(lg[lp[i]] for i in range(6))
            if pp not in group:
                group[pp]=ll;queue.append(pp)
            else:
                need(group[pp]==ll, 'Common all-root/quartic action')
    need(len(group)==46080 and len(set(group.values()))==720, 'Full actual native group')
    char=tuple(int(sum(t[5:]))%4 for t in targets)
    stabilizer=[(p,lp) for p,lp in group.items() if all(char[p[i]]==char[i] for i in range(240))]
    need(len(stabilizer)==64 and len(set(lp for p,lp in stabilizer))==16, 'Actual internal-character stabilizer')
    output={}
    for q in (2,3,4):
        pp,lp=next((p,lp) for p,lp in stabilizer if lp[5]==q and set(lp[j] for j in (0,1))=={0,1})
        G=s.Matrix.hstack(*[s.Matrix(targets[pp[ri[b]]]) for b in cs._SIMPLE])*T.inv()
        Gn=cl(O.T*G*O)
        g=s.Matrix(4,4,lambda i,j:Gn[2*i,2*j]+s.I*Gn[2*i+1,2*j])
        perm=s.zeros(6)
        for j,k in enumerate(lp):
            perm[k,j]=1
        rho=cl(N*perm*N.T)
        need(G.T*G==s.eye(8) and g.H*g==s.eye(4), 'True full source isometry')
        need(cl(rho*N-N*perm)==s.zeros(5,6), 'Literal intrinsic-ray transport')
        sub=dict(zip(x,g.T*s.Matrix(x)))
        for j,f in enumerate(polys):
            need(s.expand(f.subs(sub,simultaneous=True)-sum(rho[k,j]*polys[k] for k in range(5)))==0,
                 'True quartic representation, no fitted S6')
        need(cl(rho*polar(N,5)*rho.H-polar(N,q))==s.zeros(5), 'Actual polar mark transport')
        output[q]=(G,g,rho,lp)
    return targets,output

def main():
    ld=cs._lattice_data(); T=s.Matrix.hstack(*map(s.Matrix,cs._TARGET_SIMPLE))
    O=cl(T*ld['basis'].inv()/s.sqrt(2))
    U=s.zeros(8,4)
    for j in range(4):
        U[2*j,j]=1/s.sqrt(2);U[2*j+1,j]=-s.I/s.sqrt(2)
    M=s.Matrix(json.loads((R/'experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json').read_text())['explicit_Weyl_matrix'])
    ps,ns=old.code(); x=s.symbols('x0:4')
    polys=[s.expand(sum(old.Q.to_sympy(c)*s.prod(x[j] for j in mon) for mon,c in f.items())/s.sqrt(ns[i])) for i,f in enumerate(ps)]
    raw=s.Matrix([[2,2,-1,-1,-1,-1],[0,0,-1,-1,1,1],[0,0,-1,1,-1,1],[0,0,1,-1,-1,1],[1,-1,0,0,0,0]])
    N=cl(s.diag(2,s.sqrt(12),s.sqrt(12),s.sqrt(12),s.sqrt(24))*raw/s.sqrt(48))
    need(N*N.T==s.eye(5) and N.T*N==s.eye(6)-s.ones(6)/6, 'Original reconstructed simplex normalization')
    erasures=[s.diff(f,x[a]) for f in polys for a in range(4)]
    fields=polys+erasures
    # A finite generating bank for ALL position/mode two-point functions.
    bank=[]
    for r in range(5):
        bank.append(s.Matrix([[s.diff(f,*[x[j] for j in mon])/s.sqrt(math.prod(math.factorial(n) for n in Counter(mon).values())) if mon else f
                               for f in fields] for mon in multi(r)]))
    origin=evaluate(bank,x,s.zeros(4,1))
    need(cl((origin[3].H*origin[3])[5:,5:])==s.eye(20), 'Actual twenty erasure fields, not dimension matching')
    targets,transports=source_transports(N,polys,x,O)
    reference=json.loads((BASE/'fold_response_tensors_results.json').read_text())
    output={}
    expected={'1/4':[1200580,1300621,982285,1446340], '3/4':[1308100,1521805,1089805,1338820]}
    for alpha in (s.Rational(1,4),s.Rational(3,4)):
        delta=s.Matrix([(j+alpha)/8-s.Rational(1,2) for j in range(8)])
        d=cl(M*delta);b=cl(U.T*O.T*d)
        E=evaluate(bank,x,b)
        full=[cl(e.H*e) for e in E]
        Ds=[A[:5,:5] for A in full]
        prior=[s.Matrix(A) for A in reference['cases'][str(alpha)]['wick_tensors']]
        need(all(cl(a-bb)==s.zeros(5) for a,bb in zip(Ds,prior)), 'Independent complete original quartic tensors reproduced')
        starts={};c1s=[];plane_invariants=[]
        for q in range(2,6):
            w=s.sqrt(s.Rational(6,5))*N[:,q]
            coefficients=[s.simplify((w.T*Dr*w)[0]) for Dr in Ds]
            c1s.append(coefficients[1])
            plane_trace=s.simplify(s.trace(polar(N,q)*Ds[1]))
            plane_invariants.append(plane_trace)
            need(coefficients[1]==s.Rational(expected[str(alpha)][q-2],1006632960), 'Exact new invariant')
            need(coefficients[2:]==[s.Rational(4335,131072),s.Rational(85,256),s.Integer(1)], 'Lower-order equal controls')
            labels=[j for j in range(6) if j!=q]
            D=s.Matrix.hstack(*[N[:,j]-N[:,q] for j in labels])
            G=s.eye(5)+s.ones(5)
            invsqrt=s.eye(5)+(1/s.sqrt(6)-1)*s.ones(5)/5
            A=cl(D*invsqrt)
            need(cl(D.T*D-G)==s.zeros(5) and cl(A.T*A)==s.eye(5), 'Same original incident-field metric, exact orthonormal frame')
            need(all(cl(A[:,k]-(N[:,j]-(s.sqrt(6)-1)*N[:,q]/5))==s.zeros(5,1) for k,j in enumerate(labels)), 'Actual rooted polar columns')
            L=s.diag(A,s.kronecker_product(A,s.eye(4)))
            # Store complete mixed quartic/erasure coefficients for each q.
            joint=[cl(L.T*J*L) for J in full]
            starts[str(q)]={'labels':labels,'normalized_ray_kernel_coefficients':[str(c) for c in coefficients],
                            'polar_plane_D1_trace':str(plane_trace),
                            'rooted_orthonormal_frame':mat(A),'raw_incident_quartic_kernel':[mat(cl(D.T*J*D)) for J in Ds],
                            'full_joint_25_field_kernel_coefficients':[mat(J) for J in joint],
                            'field_degrees':[4]*5+[3]*20,
                            'mixed_kernel_rule':'Jr[a,b]/((z-w)^(2r)*z^(degree[a]-r)*w^(degree[b]-r)), sum r0..4; original i0/adjoint retained'}
        need(len(set(c1s))==4, 'All six fixed-state pairwise comparisons distinguished')
        need(all(s.simplify(a-bb)!=0 for a,bb in it.combinations(plane_invariants,2)), 'All six plane-only invariants distinct, even without distinguished root ray')
        tc={}
        for q,(G,g,rho,lp) in transports.items():
            dp=cl(G*d);bp=cl(U.T*O.T*dp)
            need(cl(bp-s.conjugate(g)*b)==s.zeros(4,1), 'Correct positive Cartan shift transport, not g*b')
            Ep=evaluate([v[:,:5] for v in bank],x,bp)
            Dp=[cl(e.H*e) for e in Ep]
            need(all(cl(rho.H*D*rho-Dr)==s.zeros(5) for D,Dr in zip(Dp,Ds)), 'All five complete kernels covariant under SAME background transport')
            z=s.Symbol('z');V=s.Matrix([f.subs(dict(zip(x,b)),simultaneous=True) for f in polys])
            Vp=s.Matrix([f.subs(dict(zip(x,bp)),simultaneous=True) for f in polys])
            need(cl(rho.T*Vp-V)==s.zeros(5,1), 'One-point phases and transport retained')
            charged=[];different=[]
            for beta2 in targets:
                beta=s.Matrix(beta2)/2;gb=G*beta
                db=s.simplify((d.T*beta)[0]);dgb=s.simplify((d.T*gb)[0]);dpair=s.simplify((dp.T*gb)[0])
                need(dpair==db, 'All 240 charged twist exponents and excitation energies covariant')
                charged.append({'root_doubled':[int(v) for v in beta2],'d_dot_root':str(db),'fixed_d_after_field_transport':str(dgb)})
                if db!=dgb:
                    different.append({'beta_doubled':[int(v) for v in beta2],'G_beta_doubled':[str(v) for v in 2*gb],
                                      'old_root_excitation':str(1+db),'fixed_state_transported_root_excitation':str(1+dgb)})
            need(bool(different), 'Actual charged source/time is not a fixed-background symmetry')
            tc[str(q)]={'logical_permutation':list(lp),'source_isometry':mat(G),'positive_native_matrix':mat(g),'quartic_matrix':mat(rho),
                        'transported_source_d':[str(v) for v in dp],'fixed_d_changed_charge_count':len(different),
                        'charged_time_witness':different[0],'all_240_charged_diagnostics':charged,
                        'covariant_full_quartic_kernel':'PASS','preserves_internal_character':True,
                        'claim':'Individual E8 lattice lift with transported cocycle cochain; no raw-CAR symmetry or globally split G31 lift assumed'}
        output[str(alpha)]={'original_d':[str(v) for v in d],'original_b':[str(v) for v in b],
                            'starts':starts,'actual_native_transports_from_q5':tc,
                            'c1_common_denominator':1006632960,'c1_numerators':expected[str(alpha)],
                            'all_pairwise_fixed_state_comparisons':'DISTINCT_WITHIN_FIXED_FIELD_REFERENCE'}
        print(json.dumps({'alpha':str(alpha),'c1_numerators':expected[str(alpha)],'six_pairs_distinct':True,
                          'complete_joint_kernel_fields':25,'source_transports':3}),flush=True)
    out={'status':'PASS_EXACT_FOUR_STARTS_STATE_RELATIVE_COMPARISON','date':'2026-10-18','cases':output,
         'all_word_completion':'Ordered E8 Heisenberg/lattice Wick generating formula in proof, including all charges, cocycles, derivatives, time and erasure fields; no executed intermediate projection',
         'verdict':'Four marked rays/flags distinguishable in SAME stipulated fold-state/current dictionary. Simultaneous complete data transport is covariance, not a fixed-state symmetry or proof of gauge redundancy. No four physical P1 theories or q selector derived.',
         'scope':'Original M/O, given E8 extension, fixed conformal NS current reference and coordinates/adjoint/time retained. Not all P1 representations; raw physical field identification, preparation and output remain open.'}
    Path(__file__).with_name('four_starting_points_results.json').write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
    print(out['status'])

def check_plane_trace_addendum():
    """Check new plane invariant against the completed FULL joint matrices.

    Separate small replay after the initial healthy full computation, whose
    already parsed version predated this strengthening. No full rerun needed.
    """
    path=Path(__file__).with_name('four_starting_points_results.json')
    out=json.loads(path.read_text())
    prior=json.loads((BASE/'fold_response_tensors_results.json').read_text())
    raw=s.Matrix([[2,2,-1,-1,-1,-1],[0,0,-1,-1,1,1],[0,0,-1,1,-1,1],[0,0,1,-1,-1,1],[1,-1,0,0,0,0]])
    N=cl(s.diag(2,s.sqrt(12),s.sqrt(12),s.sqrt(12),s.sqrt(24))*raw/s.sqrt(48))
    for alpha,case in out['cases'].items():
        D1=s.Matrix(prior['cases'][alpha]['wick_tensors'][1])
        vals=[]
        for q in range(2,6):
            item=case['starts'][str(q)]
            J1=s.Matrix(item['full_joint_25_field_kernel_coefficients'][1])
            need(item['labels'][:2]==[0,1], 'Original polar two-plane columns retained')
            actual=s.simplify(J1[0,0]+J1[1,1])
            expected=s.simplify(s.trace(polar(N,q)*D1))
            need(s.simplify(actual-expected)==0, 'Plane-only trace from full joint kernel equals independent original-D1 contraction')
            item['polar_plane_D1_trace']=str(actual)
            vals.append(actual)
        need(all(s.simplify(a-b)!=0 for a,b in it.combinations(vals,2)), 'All six exact plane-only comparisons differ')
    out['plane_only_addendum']='PASS_EXACT_ALL_12_PAIR_COMPARISONS_FROM_FULL_JOINT_KERNEL_AND_ORIGINAL_D1'
    path.write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
    print(out['plane_only_addendum'])

if __name__=='__main__':
    if '--plane-addendum' in sys.argv:
        check_plane_trace_addendum()
    else:
        main()
