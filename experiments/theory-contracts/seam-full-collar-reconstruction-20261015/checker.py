#!/usr/bin/env python3
"""Exact finite controls of a conditional, unbounded full-collar graph theorem."""
from __future__ import annotations
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sympy as s

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[2]


def canonical(value):
    return json.dumps(value,sort_keys=True,separators=(',',':'))


def hash_value(value):
    return sha256(canonical(value).encode()).hexdigest()


def exact(matrix):
    return [[str(item) for item in row] for row in matrix.tolist()]


def zero(matrix):
    return matrix==s.zeros(*matrix.shape)


def power_two_exponent(value):
    value=Fraction(value)
    numerator,denominator=value.numerator,value.denominator
    if numerator<=0 or numerator&(numerator-1) or denominator&(denominator-1):
        raise ValueError('Diagnostic inversion requires an exact positive power of two')
    return numerator.bit_length()-denominator.bit_length()


def rational_entry(item):
    if not isinstance(item,(str,int)) or isinstance(item,bool):
        raise ValueError('Matrix entries must be exact rational strings or integers')
    q=Fraction(item)
    return s.Rational(q.numerator,q.denominator)


def graph(transfer):
    identity=s.eye(transfer.rows)
    a=(identity+transfer.T*transfer).inv()
    c=transfer*a
    return a,c,a.row_join(c).col_join(c.row_join(identity-a))


def reconstruct_kernel(document):
    """Consume only the supplied two-face graph, never a supplied B or spectrum."""
    frame='ordered_shared_orthonormal_frame'
    if document.get('fiber_identification')!=frame:
        raise ValueError('Explicit fiber_identification=ordered_shared_orthonormal_frame is required: positivity and selfadjointness of transfer need the given common orthonormal fiber comparison')
    faces=document['boundary_faces']
    if len(faces)!=2 or faces[0]['label']==faces[1]['label']:
        raise ValueError('Two distinct labelled boundary faces are required')
    n=len(faces[0]['modes'])
    if not n or len(faces[1]['modes'])!=n:
        raise ValueError('Both explicitly labelled faces must have the same positive dimension')
    if any(len({canonical(mode) for mode in face['modes']})!=n for face in faces):
        raise ValueError('Mode labels on each face must be distinct')
    raw=document['graph_projection']
    if len(raw)!=2*n or any(len(row)!=2*n for row in raw):
        raise ValueError('Graph dimension disagrees with the two boundary mode lists')
    g=s.Matrix([[rational_entry(item) for item in row] for row in raw])
    if g.T!=g or g*g!=g:
        raise ValueError('Supplied graph matrix must be an exact orthogonal projection')
    a=g[:n,:n]
    if a.det()==0:
        raise ValueError('C11 is singular; sign-only data does not determine a finite collar transfer')
    t=g[n:,:n]*a.inv()
    if t.det()==0:
        raise ValueError('Transfer is not invertible; sign-only data is not zero collar time')
    if t.T!=t or any(t[:k,:k].det()<=0 for k in range(1,n+1)):
        raise ValueError('Reconstructed transfer must be selfadjoint and strictly positive')
    if graph(t)[2]!=g:
        raise ValueError('Supplied projection has incorrect product-collar graph blocks')
    length=document.get('normal_length')
    if isinstance(length,dict):
        if length.get('basis')!='log(2)':
            raise ValueError('Supported exact length bases are rational and coefficient times log(2)')
        coefficient=rational_entry(length['coefficient'])
        ell=coefficient*s.log(2)
    else:
        coefficient=rational_entry(length)
        ell=coefficient
    if coefficient<=0:
        raise ValueError('Given normal collar length must be strictly positive')
    return {'outcome':'CONDITIONAL_OPERATOR_RECONSTRUCTED_FROM_SUPPLIED_GRAPH',
        'dimension':n,'boundary_faces':faces,'fiber_identification':frame,'normal_length':str(ell),
        'transfer':exact(t),'logit_argument':exact(a*(s.eye(n)-a).inv()),
        'occupation_from_second_face':exact(g[n:,n:]),
        'state_time_scope':'If the field is charged CAR on both faces and has covariance G, this n=G22 is derived from that same datum, not independently fitted; log((I-n)/n)=2 ell B.',
        'operator':{'operation':'negative_selfadjoint_spectral_log',
                    'argument':exact(t),'divide_by_given_length':str(ell)},
        'operator_meaning':'B=-log(T)/ell, uniquely defined by positive spectral functional calculus; no B was supplied to the reader',
        'source_scope':'conditional exact finite product-collar graph on given labelled faces and length; no identification with an actual P1 matrix kernel',
        'physical_source_selection_proved':False}


def audit_carrier_kernel(document):
    """Validate an actual supplied carrier embedding, never invent a P2 split."""
    base=reconstruct_kernel(document)
    g=s.Matrix([[rational_entry(item) for item in row] for row in document['graph_projection']])
    raw_w=document['carrier_embedding']
    raw_e=document['claimed_carrier_involution']
    k=len(raw_w[0]) if raw_w else 0
    if not k or len(raw_w)!=g.rows or any(len(row)!=k for row in raw_w):
        raise ValueError('Carrier embedding must have 2n rows and a positive, consistent carrier dimension')
    if len(raw_e)!=k or any(len(row)!=k for row in raw_e):
        raise ValueError('Claimed carrier involution dimension disagrees with the embedding')
    w=s.Matrix([[rational_entry(item) for item in row] for row in raw_w])
    e=s.Matrix([[rational_entry(item) for item in row] for row in raw_e])
    identity=s.eye(k)
    if w.T*w!=identity:
        raise ValueError('Carrier embedding is not isometric: W*W must equal I')
    if e.T!=e or e*e!=identity:
        raise ValueError('Claimed carrier operator must be selfadjoint and involutive')
    f=2*g-s.eye(g.rows)
    positive_only=g*w==w
    if positive_only and e!=identity:
        raise ValueError('Carrier lies in Ran(G): (2G-I)W=W, so injective W forces E=I. A negative carrier label after positive-only compression is incompatible with the supplied source dictionary.')
    if f*w!=w*e:
        raise ValueError('Carrier dictionary fails actual intertwining: (2G-I)W is not W E')
    return {**base,'outcome':'GIVEN_CARRIER_DICTIONARY_VERIFIED_NOT_DERIVED',
        'carrier_audit':{'dimension':k,'embedding_isometric':True,
            'claimed_involution_selfadjoint_and_involutive':True,
            'actual_intertwining_verified':True,'image_in_Ran_G':positive_only,
            'actual_compression_W_star_F_W':exact(w.T*f*w),
            'both_polarizations_actually_carried':e!=identity and e!=-identity,
            'positive_only_implication':'If G W=W then (2G-I)W=W and W*(2G-I)W=I; an intertwining E must be I.',
            'scope':'The supplied embedding and dictionary were checked, not derived from P1/P2. A carrier spanning both polarizations can support a conditional internal marking; no 3+2 carrier or low-energy five-mode choice is manufactured.'}}


def run_checks():
    modes=list(range(-2,2))
    quarter=s.Rational(1,4)
    h=s.diag(*(n+quarter for n in modes))
    identity=s.eye(h.rows)
    null=s.zeros(h.rows)
    b=s.diag(h,-h)
    j=null.row_join(-identity).col_join(identity.row_join(null))
    reflection=null.row_join(identity).col_join(identity.row_join(null))
    size=b.rows
    t=s.diag(*(s.Rational(2)**(-4*b[k,k]) for k in range(size)))
    a,c,g=graph(t)
    full_identity=s.eye(2*size)
    reconstructed_t=c*a.inv()
    reconstructed_b=s.diag(*(s.Rational(power_two_exponent(a[k,k]/(1-a[k,k])),8)
                             for k in range(size)))
    occupation=s.eye(size)-a
    occupation_ratio=(s.eye(size)-occupation)*occupation.inv()
    occupation_b=s.diag(*(s.Rational(power_two_exponent(occupation_ratio[k,k]),8)
                         for k in range(size)))
    half=s.diag(*(1 if b[k,k]>0 else 0 for k in range(size)))
    boundary_reflection=s.zeros(size).row_join(reflection).col_join(
        reflection.row_join(s.zeros(size)))
    compact_deck=s.diag(reflection,reflection)
    checks={
        'graph_selfadjoint':g.T==g,
        'graph_idempotent':g*g==g,
        'graph_range_relation':zero(g[size:,:]-t*g[:size,:]),
        'graph_orthogonal_complement':zero(g*(-t).col_join(s.eye(size))),
        'transfer_reconstructed_exactly':reconstructed_t==t,
        'B_reconstructed_exactly_without_floating_log':reconstructed_b==b,
        'reflection_involution':reflection*reflection==s.eye(size),
        'Clifford_J_square_minus_one':j*j==-s.eye(size),
        'reflection_reverses_J':reflection*j*reflection==-j,
        'reflection_reverses_B':reflection*b*reflection==-b,
        'J_anticommutes_B':zero(j*b+b*j),
        'reflection_exchanges_half_polarizations':reflection*half*reflection==s.eye(size)-half,
        'reflection_does_not_preserve_half_polarization':reflection*half*reflection!=half,
        'full_two_face_reflection_preserves_graph':boundary_reflection*g*boundary_reflection==g,
        'compact_deck_fixes_each_boundary_but_not_graph':compact_deck*g*compact_deck!=g,
        'compact_deck_exchanges_diagonal_graph_blocks':
            compact_deck*g*compact_deck==(s.eye(size)-a).row_join(c).col_join(c.row_join(a)),
        'compact_deck_does_not_equal_graph_polarization':compact_deck!=2*g-full_identity,
        'occupation_derived_from_same_graph':occupation==g[size:,size:],
        'occupation_eigenvalues_strictly_between_zero_and_one':
            all(0<occupation[k,k]<1 for k in range(size)),
        'occupation_faithful_modular_ratio':occupation_ratio==t.inv()**2,
        'same_B_recovered_from_occupation':occupation_b==b,
    }
    # Mix opposite signs: a nontrivial frame change must transport the state too.
    u=s.eye(size)
    u[0,0]=u[2,2]=s.Rational(3,5)
    u[0,2]=-s.Rational(4,5)
    u[2,0]=s.Rational(4,5)
    framed_t=u*t*u.T
    framed_a,framed_c,framed_g=graph(framed_t)
    doubled_u=s.diag(u,u)
    marker=s.diag(*range(1,size+1))
    marker[0,1]=marker[1,0]=1
    framed_b=u*b*u.T
    framed_marker=u*marker*u.T
    framed_occupation=s.eye(size)-framed_a
    checks.update({
        'frame_unitary':u*u.T==s.eye(size),
        'frame_changes_polarization':u*half*u.T!=half,
        'frame_B_is_not_diagonal':framed_b!=s.diag(*framed_b.diagonal()),
        'full_graph_frame_covariance':framed_g==doubled_u*g*doubled_u.T,
        'A_frame_covariance':framed_a==u*a*u.T,
        'transfer_frame_reconstruction':framed_c*framed_a.inv()==framed_t,
        'B_logit_spectral_calculus_frame_covariance':
            sum((s.Rational(power_two_exponent(a[k,k]/(1-a[k,k])),8)
                 *(u[:,k]*u[:,k].T) for k in range(size)),s.zeros(size))==framed_b,
        'marked_source_commutator_covariance':
            framed_b*framed_marker-framed_marker*framed_b==u*(b*marker-marker*b)*u.T,
        'marker_commutator_is_nonzero':not zero(b*marker-marker*b),
        'state_frame_covariance':s.eye(size)-u*half*u.T==u*(s.eye(size)-half)*u.T,
        'faithful_occupation_frame_covariance':framed_occupation==u*occupation*u.T,
        'framed_occupation_positive_and_co_positive':all(
            framed_occupation[:k,:k].det()>0 and (s.eye(size)-framed_occupation)[:k,:k].det()>0
            for k in range(1,size+1)),
        'occupation_modular_ratio_frame_covariance':
            (s.eye(size)-framed_occupation)*framed_occupation.inv()==framed_t.inv()**2,
        'same_framed_B_from_occupation_spectral_log':
            sum((s.Rational(power_two_exponent(occupation_ratio[k,k]),8)
                 *(u[:,k]*u[:,k].T) for k in range(size)),s.zeros(size))==framed_b,
    })
    labels=[]
    for n in modes:
        labels.append({'n':n,'alpha':'1/4','negative_partner_n':-n-1,'negative_partner_alpha':'3/4',
                       'b':str(n+quarter),'negative_b':str(-n-quarter)})
    checks['opposite_quarter_domains_preserved']=all(
        s.Rational(row['negative_partner_n'])+s.Rational(3,4)==-s.Rational(row['b']) for row in labels)
    # Reader gets a graph matrix only: T/B are not input fields.
    fixture={'fiber_identification':'ordered_shared_orthonormal_frame',
             'boundary_faces':[{'label':'r=0','modes':['a','b']},
                               {'label':'r=ell','modes':['a','b']}],
             'normal_length':'1','graph_projection':exact(graph(s.Matrix([[2,1],[1,2]]))[2])}
    recovered=reconstruct_kernel(fixture)
    checks['generic_reader_reconstructs_from_graph_only']=recovered['transfer']==[['2','1'],['1','2']]
    rejections={}
    for label,invalid in {
        'sign_only':s.diag(1,1,0,0),
        'opposite_sign_only':s.diag(0,0,1,1),
        'negative_transfer':graph(-s.eye(2))[2],
        'not_projection':s.eye(4)/2,
    }.items():
        try:
            reconstruct_kernel({**fixture,'graph_projection':exact(invalid)})
        except ValueError as error:
            rejections[label]=str(error)
        else:
            raise ArithmeticError('Invalid kernel accepted: '+label)
    checks['reader_rejects_sign_limits_and_nonpositive_or_invalid_graphs']=len(rejections)==4
    frame_rejections={}
    for label,value in (('missing',None),('different','unrelated_boundary_frames')):
        invalid={**fixture}
        if value is None: invalid.pop('fiber_identification')
        else: invalid['fiber_identification']=value
        try: reconstruct_kernel(invalid)
        except ValueError as error: frame_rejections[label]=str(error)
        else: raise ArithmeticError('Missing or different fiber identification accepted')
    checks['reader_rejects_missing_or_different_fiber_identification']=len(frame_rejections)==2
    carrier_fixture={
        'fiber_identification':'ordered_shared_orthonormal_frame',
        'boundary_faces':[{'label':'r=0','modes':['synthetic']},{'label':'r=ell','modes':['synthetic']}],
        'normal_length':'1','graph_projection':[['16/25','12/25'],['12/25','9/25']],
        'carrier_embedding':[['4/5'],['3/5']],'claimed_carrier_involution':[['1']]}
    positive_audit=audit_carrier_kernel(carrier_fixture)
    checks['carrier_positive_only_identity_dictionary_accepted']=positive_audit['carrier_audit']['image_in_Ran_G']
    mixed_fixture={**carrier_fixture,'carrier_embedding':[['4/5','-3/5'],['3/5','4/5']],
                   'claimed_carrier_involution':[['1','0'],['0','-1']]}
    mixed_audit=audit_carrier_kernel(mixed_fixture)
    checks['carrier_both_polarizations_actual_intertwiner_accepted']=mixed_audit['carrier_audit']['both_polarizations_actually_carried']
    carrier_rejections={}
    for label,mutation in {
        'negative_label_on_positive_only_carrier':{'claimed_carrier_involution':[['-1']]},
        'nonisometric_embedding':{'carrier_embedding':[['4/5'],['0']]},
        'noninvolutive_label':{'claimed_carrier_involution':[['1/2']]},
        'nonintertwining_embedding':{'carrier_embedding':[['1'],['0']]},
    }.items():
        try: audit_carrier_kernel({**carrier_fixture,**mutation})
        except ValueError as error: carrier_rejections[label]=str(error)
        else: raise ArithmeticError('Invalid carrier dictionary accepted: '+label)
    checks['carrier_reader_rejects_invalid_dictionaries']=len(carrier_rejections)==4
    checks['carrier_negative_label_has_precise_positive_compression_reason']='forces E=I' in carrier_rejections['negative_label_on_positive_only_carrier']
    if not all(checks.values()):
        raise ArithmeticError('Exact control failed: '+str([k for k,v in checks.items() if not v]))
    return {
        'outcome':'EXACT_CONDITIONAL_FULL_COLLAR_RECONSTRUCTION',
        'finite_diagnostic':{'modes_per_branch':modes,'branches':['1/4','3/4'],
            'length':'4 log(2)','B_dimension':size,'graph_dimension':2*size,
            'rational_transfer_exponents':[-int(4*b[k,k]) for k in range(size)],
            'quarter_labels':labels},
        'checks':checks,'check_count':len(checks),
        'reader_fixture':fixture,'reader_output':recovered,'reader_rejection_controls':rejections,
        'fiber_identification_rejection_controls':frame_rejections,
        'carrier_positive_fixture':carrier_fixture,'carrier_positive_audit':positive_audit,
        'carrier_mixed_audit':mixed_audit,'carrier_rejection_controls':carrier_rejections,
        'exact_data':{'B':exact(b),'T':exact(t),'A':exact(a),'C21':exact(c),
                      'occupation_from_G22':exact(occupation),'occupation_modular_ratio':exact(occupation_ratio),
                      'occupation_eigenvalues':[str(occupation[k,k]) for k in range(size)],
                      'framed_occupation':exact(framed_occupation),
                      'frame_U':exact(u),'framed_B':exact(framed_b)},
        'analytic_statement':{
            'premises':'B selfadjoint, known length ell>0, exact normal-independent product equation (partial_r+B)f=0, complete L2 two-face Cauchy data',
            'T':'exp(-ell B), positive injective densely defined closed operator; D(T)={u: integral exp(-2 ell lambda) dmu_u<infinity}',
            'Cauchy_space':'closed graph {(u,T u): u in D(T)} in H direct_sum H',
            'graph_projection':'[[A,C],[C,1-A]], A=(1+exp(-2 ell B))^-1, C=(2 cosh(ell B))^-1',
            'B_inverse':'B=(2 ell)^-1 logit(A), via spectral functional calculus',
            'B_domain':'{u: integral |log(a/(1-a))|^2 dmu_A,u(a)<infinity}',
            'T_inverse':'T=closure(C A^-1), or joint functional calculus C/A on D(T)',
            'raw_composition_domain':'D(C A^-1)=Ran(A)=D(T^2), generally strictly smaller than D(T)',
            'strict_domain_example':'B e_n=-n e_n, n>=1; u_n=exp(-3 ell n/2) belongs to D(T) but not D(T^2)',
            'endpoint_spectrum':'A and 1-A injective with dense range; 0 and 1 may be continuous spectrum, never eigenvalues',
            'length_unknown':'Only ell B is recoverable if ell is not independently given',
            'trace_scope':'L2 endpoint Cauchy space is the graph closure of smooth spectral-core solutions, not an assertion that every trace gives an H1 solution up to both faces',
            'reflection':'R f(ell-r), RJR=-J and RBR=-B, commutes with J(partial_r+B); on traces (u,v)->(Rv,Ru)',
            'compact_double_counterexample':'Double [0,ell] x S1 to an r-circle of length 2ell with the same constant product Dirac operator. The deck r->-r fixes both boundary faces and acts on traces as diag(R,R), which exchanges A with I-A and does not preserve G. It is distinct from the inside-cylinder reflection swapping the two faces.',
            'faithful_state_time_link':'Given charged CAR(H direct_sum H) and pure quasifree covariance G, restriction to the second face has n=G22=I-A=(1+exp(2ell B))^-1. Its modular ratio (I-n)/n=T^-2 recovers the same B, with no independently supplied occupation matrix.',
            'state_scope':'CAR and the covariance-state physical identification remain additional premises. The infinite-dimensional KMS state is not asserted to be a trace-class Gibbs density matrix in an arbitrary representation.',
            'quarter_frame_scope':'R=sigma_x acts on the common reference Hilbert space for diag(A,-A). Opposite spectral shifts are relabelled -(n+1/4)=(-n-1)+3/4; literal twisted function domains need an explicit pullback identification, not a scalar identity map.',
            'half_collar_limit':'A tends strongly to 1_(0,infinity)(B)+1/2*1_{0}(B) as ell grows; in the nonzero-quarter diagnostic it loses all spectral magnitudes',
            'frame_scope':'Finite exact simultaneous transport of graph, B, marker and state; an arbitrary frame rotation is not claimed pointwise local in a fixed geometric chart'},
        'unproved_original_source_obligations':[
            'P1 identifies the actual selfadjoint local boundary/collar operator and common field space',
            'original kernel supplies the complete finite two-face Cauchy graph, not just its sign limit',
            'original source supplies or calibrates ell in the same normal scale',
            'exact product form holds for the whole retained collar, not merely its principal symbol',
            'physical CAR, reference state, P2 carrier dictionary and events are identified'],
        'physical_source_selection_proved':False,
    }


def package():
    paths=[Path(__file__).resolve(),REPO/'tfpt_1_architecture_e8.tex',
           REPO/'_archive/tfpt-45/source_extracts/01_boundary_kernel_source.tex',
           REPO/'_archive/tfpt-45/source_extracts/02_carrier_source.tex']
    for name in ('input-example.json','PROOF.txt','atlas_snapshot.json','lean/FullCollar.lean','lean-receipt.json'):
        if (HERE/name).exists(): paths.append(HERE/name)
    payload=run_checks()
    return {'schema':'full-collar-graph-reconstruction-v1','payload':payload,
            'payload_sha256':hash_value(payload),
            'sources':{path.relative_to(REPO).as_posix():sha256(path.read_bytes()).hexdigest() for path in paths}}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=('check','verify','reconstruct','audit-carrier'),nargs='?',default='check')
    parser.add_argument('--kernel',type=Path)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    if args.output is None:
        args.output=HERE/({'reconstruct':'reconstruction.json','audit-carrier':'carrier-audit.json'}.get(args.command,'results.json'))
    if args.command in ('reconstruct','audit-carrier'):
        if args.kernel is None: parser.error(args.command+' requires --kernel input.json')
        reader=audit_carrier_kernel if args.command=='audit-carrier' else reconstruct_kernel
        result=reader(json.loads(args.kernel.read_text()))
        result['input_file_sha256']=sha256(args.kernel.read_bytes()).hexdigest()
        result['input_path']=str(args.kernel.resolve())
        args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
        print(json.dumps(result,sort_keys=True))
        return
    result=package()
    if args.command=='check':
        args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    elif json.loads(args.output.read_text())!=result:
        raise ValueError('Exact replay/source mismatch')
    print(json.dumps({'outcome':result['payload']['outcome'],'checks':result['payload']['check_count'],
                      'payload_sha256':result['payload_sha256'],'verified':args.command=='verify'}))


if __name__=='__main__': main()
