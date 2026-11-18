"""Original polar mark, erasure decoder and the same marked-clock formula."""
from pathlib import Path
import importlib.util
import json
import sys
import sympy as s

REPO = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
BASE = REPO / 'experiments/theory-contracts/joint-native-algebra-selection-20261016'
sys.path.insert(0, str(REPO))
spec = importlib.util.spec_from_file_location('original_native_audit', BASE / 'operator_checker.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)
from tfpt_explorer.marker_selection import build_marker_data
from tfpt_explorer.charge_selection import _parse_matrix, _process_frame
from tfpt_explorer.marked_source_process import _source_operators
from tfpt_explorer.process import _code_basis, _invariants
from tfpt_explorer.cartan_source import _lattice_data, _rays

clean = audit.clean
encode = audit.matrix_record

def require(condition, label):
    if not condition:
        raise ArithmeticError(label)

def main():
    data = build_marker_data()['data']
    _, supports = _code_basis()
    native_sigma, _, _ = _source_operators()
    sigma5 = audit.quartic_compression(native_sigma, supports)
    rotation = _parse_matrix(data['basis_alignment']['rotation_matrix'])
    process_frame = _process_frame()
    # The original W6 simplex is N=O C^T; its labels reverse the C labels.
    # Retain the ORIGINAL native four-clock and use N consistently throughout.
    frame = clean(rotation * process_frame.T).T
    def logical(p):
        return clean(frame.T * audit.permutation_matrix(tuple(p)) * frame)
    sigma_expected = logical((0, 1, 3, 4, 2, 5))
    require(sigma5 == sigma_expected, 'Actual native sigma must match the original quotient')
    reverse = audit.permutation_matrix(tuple(reversed(range(6))))
    require(clean(rotation + process_frame.T*reverse*process_frame)==s.zeros(5), 'O is minus the native reverse permutation')
    y5 = _parse_matrix(data['native_polar_transport']['native_charge_matrix'])
    p2 = clean((6*y5 + 2*s.eye(5))/5)
    p3 = s.eye(5)-p2
    tau5 = logical((1, 0, 2, 3, 4, 5))
    f5 = clean((s.eye(5)+sigma5+sigma5**2)/3)
    r5 = clean(p3*(s.eye(5)-2*f5)+s.I*p2*tau5)
    require(clean(r5.H*r5)==s.eye(5), 'Original R5 unitarity')

    # The original decoder, not an arbitrary twenty-dimensional identification.
    w = audit.erasure_isometry(supports)
    pi = clean(w.H*w)
    require(pi == s.eye(20), 'Actual erasure isometry')
    y20 = s.kronecker_product(y5, s.eye(4))
    p220 = s.kronecker_product(p2, s.eye(4))
    y3 = clean(w*y20*w.H)
    require(clean(w.H*y3*w)==y20, 'Original three-register decoder action')
    r20_decoded = s.kronecker_product(r5, s.eye(4))
    # The unit of the decoder corner is WW*, with the complement left unchanged.
    r3 = clean(w*r20_decoded*w.H+s.eye(64)-w*w.H)
    require(clean(r3.H*r3)==s.eye(64), 'Decoded three-register unitary')

    inv = _invariants()
    roots = _rays(_lattice_data()['roots'])
    # Source N labels (0,1) are actual process C labels (5,4).
    tau_pair_id = inv['transposition_pairs'].index((4, 5))
    ids = [i for i,l in enumerate(inv['labels']) if inv['logical_to_pair'][l]==tau_pair_id]
    require(len(ids)==4, 'All four actual primitive lifts of the weak transposition')
    sigma20 = s.kronecker_product(sigma5, native_sigma.conjugate())
    f20 = clean((s.eye(20)+sigma20+sigma20**2)/3)
    s23 = logical((0, 1, 3, 2, 4, 5))
    rho = clean((s.eye(5)-f5)*(s.eye(5)+s23)/2)
    require(rho.rank()==1 and s.trace(rho)==1 and rho*rho==rho, 'Actual pure color-standard witness')
    erased = s.kronecker_product(rho,s.eye(4)/4)
    target = s.kronecker_product(clean(r5*rho*r5.H),s.eye(4)/4)
    require(clean(target-erased)==s.zeros(20), 'Logical R5 fixes the witness')
    rows=[]
    for index in ids:
        tau4=audit.native_reflection(roots[index])
        require(audit.quartic_compression(tau4,supports)==tau5,'Same actual weak-transposition quotient')
        tau20=s.kronecker_product(tau5,tau4.conjugate())
        same_polynomial=clean((s.eye(20)-p220)*(s.eye(20)-2*f20)+s.I*p220*tau20)
        require(clean(same_polynomial.H*same_polynomial)==s.eye(20),'Same-polynomial20 unitarity')
        require(clean(same_polynomial**2)==s.eye(20)-2*p220,'Same-polynomial20 square')
        delta=clean(same_polynomial*erased*same_polynomial.H-target)
        hs2=s.simplify(s.trace(delta.H*delta))
        eig=delta.eigenvals()
        distance=s.simplify(sum(abs(v)*m for v,m in eig.items())/2)
        require(hs2==s.Rational(1,4) and distance==s.Rational(1,2),'Exact common-erasure mismatch')
        rows.append({'native_ray_id':index,'unitary':True,'square_correct':True,'erasure_difference_rank':delta.rank(),'erasure_difference_HS_squared':str(hs2),'erasure_trace_distance':str(distance),'eigenvalues':{str(k):v for k,v in eig.items()}})

    # Mark q* is a retained input. Test all original sigma-fixed choices, using
    # the same positive polar prescription, not a fit to the desired spectrum.
    raw=s.Matrix(data['native_polar_transport']['raw_simplex_W6'])
    metric=s.diag(*data['native_polar_transport']['raw_metric'])
    sqrt_metric=s.diag(*[s.sqrt(v) for v in data['native_polar_transport']['raw_metric']])
    polar_factor=s.eye(5)+(1/s.sqrt(6)-1)*s.ones(5)/5
    marks={}
    for q in (0,1,5):
        others=[i for i in range(6) if i!=q]
        diff=raw[:,others]-raw[:,q]*s.ones(1,5)
        transport=clean(sqrt_metric*diff*polar_factor/s.sqrt(48))
        require(clean(transport.H*transport)==s.eye(5),'Original positive polar isometry')
        slot=s.diag(*[1 if i in (0,1,5) else 0 for i in others])
        mark=clean(transport*slot*transport.H)
        require(mark.rank()==2 and clean(mark*mark-mark)==s.zeros(5),'Conditional3+2 from sigma slot partition')
        marks[q]=mark
    require(clean(marks[5]-p2)==s.zeros(5),'Reconstruct original P2 with original polar prescription')
    simplex_native = clean(sqrt_metric*raw/s.sqrt(48))
    moving_barycenter = sum((simplex_native[:,j] for j in (2,3,4)),s.zeros(5,1))
    polar_color_lines={}
    for q in marks:
        color_line=clean((moving_barycenter-3*(s.sqrt(6)-1)*simplex_native[:,q]/5)/s.sqrt(3))
        require(clean(color_line.H*color_line)==s.ones(1,1),'Derived polar color line is normalized')
        require(clean(f5-color_line*color_line.H-marks[q])==s.zeros(5),'Source formula P2=F_sigma-c_q c_q*')
        polar_color_lines[q]=encode(color_line)

    centralizing_lifts=[]
    for pair in ((0,1),(0,5),(1,5)):
        process_pair=tuple(sorted((5-pair[0],5-pair[1])))
        pid=inv['transposition_pairs'].index(process_pair)
        candidates=[i for i,l in enumerate(inv['labels']) if inv['logical_to_pair'][l]==pid]
        p=list(range(6));p[pair[0]],p[pair[1]]=p[pair[1]],p[pair[0]]
        logical_swap=logical(p)
        commuting=[]
        for i in candidates:
            g4=audit.native_reflection(roots[i])
            require(audit.quartic_compression(g4,supports)==logical_swap,'Exact fixed-point transposition lift')
            if clean(g4*native_sigma-native_sigma*g4)==s.zeros(4):commuting.append(i)
        require(len(commuting)==1,'Unique primitive clock-centralizing lift in each four-fibre')
        equivariant=all(clean(logical_swap*marks[q]*logical_swap.H-marks[p[q]])==s.zeros(5) for q in marks)
        require(equivariant,'Original mark family covariant under actual compiler-clock centralizer')
        centralizing_lifts.append({'fixed_point_pair':list(pair),'native_reflection_ids':commuting,'polar_mark_family_equivariant':True})

    result={'status':'PASS_EXACT_MARKED_DECODER_AND_COMMON_POLYNOMIAL_TEST','frame_transport':{'N_equals_O_CT':True,'O_equals_minus_rho_reverse':True,'native_sigma_in_N_labels':[0,1,3,4,2,5],'native_sigma_in_C_labels':[0,3,1,2,4,5],'warning':'Transporting Y from N to C while retaining N-labelled permutations changes the native clock unless its four-register lift is transported too. This checker retains original N labels and the unchanged original native four-clock.'},'original_mark':{'P2_reconstructed_exactly':True,'source_fixed_point_choices':[0,1,5],'all_rank2':True,'polar_color_lines':polar_color_lines,'derived_formula':'P2(q)=F_sigma-c_q c_q*, c_q=(sum_{j=2,3,4} n_j - 3(sqrt(6)-1)n_q/5)/sqrt(3)','native_clock_centralizing_lifts':centralizing_lifts,'choice_scope':'One orbit under the native compiler-clock centralizer. This is not a proof that every anisotropic source process or physical state has this symmetry.'},'actual_decoder':{'W_shape':list(w.shape),'Y3_from_original_erasure_decoder':True,'R3_unitary_extension':True,'R20_action':'R5 tensor I4, derived as encoded logical decoder action','scope':'The decoder specifies a logical operation on the retained code; it does not identify that operation with a primitive native intervention.'},'same_polynomial_native_lift':{'sigma20':'sigma5 tensor conjugate(sigma4)','tau20':'tau5 tensor conjugate(tau4), all four primitive lifts tested','P220':'P2 tensor I4 from original erasure decoder','rows':rows,'conclusion':'Same polynomial is unitary and has the correct square but fails the full branchwise erasure action on every actual weak-transposition lift.'},'witness':{'rho5':encode(rho),'trace_distance_all_four':'1/2'},'not_claimed':['No physical primitive intervention rule or unique source state derived.','No contradiction in the existing logical clock.','No exclusion of collective logical operations or every common algebra completion.']}
    out=Path(__file__).with_name('marked_lift_results.json');out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':result['status'],'clock_lift_trace_distance':[x['erasure_trace_distance'] for x in rows],'compiler_centralizer':centralizing_lifts},indent=2))

if __name__=='__main__':main()
