"""Exact finite-source checks. Supplied original operations remain supplied inputs.

No physical source identity, regulator removal or measure is inferred from a
role label. PASS concerns the represented exact inputs and the stated branch.
Missing original operations remain MISSING_INPUT. No per-field EOM filter.
"""
import sympy as s
from pathlib import Path
import hashlib,json

MARKS=('Z0','ZU','ZY','ZUY')
IDS=('original_joint_weighted_Z0_ZU_ZY_ZUY','original_OS_history_RP',
     'actual_quantum_Ward_SD_contacts','whole_physical_P_Ref_Eta_domains',
     'original_regulator_limit_order')
class MissingInput(Exception):pass

def need(d,key):
    if key not in d or d[key] is None:raise MissingInput(key)
    return d[key]
def scalar(x):
    if isinstance(x,(list,tuple)) and len(x)==2:return s.Rational(str(x[0]))+s.I*s.Rational(str(x[1]))
    return s.sympify(x)
def matrix(x):
    if isinstance(x,s.MatrixBase):return x
    return s.Matrix([[scalar(z)for z in row]for row in x])
def clean(x):
    return x.applyfunc(s.simplify) if isinstance(x,s.MatrixBase) else s.simplify(x)
def zero(x):
    vals=list(x) if isinstance(x,s.MatrixBase) else [x]
    states=[s.simplify(v).is_zero for v in vals]
    return False if False in states else (True if all(v is True for v in states) else None)
def equal(a,b):return zero(clean(a-b))
def nonnegative(x):return s.simplify(x).is_nonnegative
def psd(A):
    hermitian=equal(A,A.H)
    if hermitian is not True:return hermitian
    # Hermitian roots are real. All elementary symmetric coefficients >=0
    # exclude every negative root via prod(t+lambda_i)>0 for t>0.
    coefficients=A.charpoly().all_coeffs()[1:]
    signs=[nonnegative((-1)**(k+1)*v)for k,v in enumerate(coefficients)]
    return False if False in signs else (True if all(v is True for v in signs) else None)
def row(cid,tests,**extra):
    status='FAIL'if False in tests.values() else ('PASS'if all(v is True for v in tests.values())else'UNDECIDED_EXACT_SIGN')
    return {'id':cid,'status':status,'tests':tests,'scope':'exact represented finite data; original-source identification is not inferred',**extra}
def execute(cid,fn,d):
    try:return fn(d)
    except MissingInput as e:return {'id':cid,'status':'MISSING_INPUT','missing_inputs':[str(e)],'physical_source_complete':None}
    except (ValueError,TypeError,KeyError,s.ShapeError)as e:return {'id':cid,'status':'INVALID_INPUT','error':str(e)}
def jsonable(x):
    if isinstance(x,s.MatrixBase):return [[str(v)for v in list(x.row(i))]for i in range(x.rows)]
    if isinstance(x,dict):return {str(k):jsonable(v)for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [jsonable(v)for v in x]
    if isinstance(x,s.Basic):return str(x)
    return x

def quotient(d,forms):
    G=matrix(need(d,'gram_form'));V=matrix(need(d,'quotient_lift'));N=matrix(need(d,'null_basis'))
    if not G.rows or not V.cols:raise MissingInput('nonempty original retained quotient')
    tests={'gram_positive':psd(G),'quotient_metric_identity':equal(V.H*G*V,s.eye(V.cols)),
           'declared_null_is_null':zero(G*N),'quotient_and_null_span_whole_basis':V.row_join(N).rank()==G.rows}
    converted={}
    for name,A in forms.items():
        A=matrix(A);tests[name+'_left_null_descent']=zero(N.H*A);tests[name+'_right_null_descent']=zero(A*N)
        converted[name]=clean(V.H*A*V)
    return tests,converted

def common_four(d):
    def impl(d):
        atoms=need(d,'atoms');source={key:None for key in MARKS};identifiers=[]
        if not atoms:raise MissingInput('nonempty original atom list')
        for atom in atoms:
            identifiers.append(need(atom,'original_configuration_id'))
            m=scalar(need(atom,'original_measure_mass'));w=scalar(need(atom,'unnormalized_scalar_weight'))
            forms=need(atom,'ordered_source_forms')
            for k in MARKS:
                A=m*w*matrix(need(forms,k));source[k]=A if source[k]is None else source[k]+A
        tests,Q=quotient(d,{**source,'P':need(d,'physical_selector_form'),'cap':need(d,'boundary_pairing_form')})
        P=Q['P'];cap=Q['cap'];ref=scalar(need(d,'reference_scalar'));outside=scalar(need(d,'scalar_factor'))
        tests.update(unique_original_atoms=len(set(identifiers))==len(identifiers),selector_adjoint=equal(P,P.H),selector_idempotent=equal(P*P,P),nonzero_absolute_reference=(ref.is_zero is False))
        if ref.is_zero is not False:return row(IDS[0],tests,values={k:None for k in MARKS})
        values={k:clean(outside*s.trace(cap*P*Q[k]*P)/ref)for k in MARKS}
        if d.get('observed_values')is not None:
            for k in MARKS:tests[k+'_matches_observed']=equal(values[k],scalar(need(d['observed_values'],k)))
        normalized=None
        if values['Z0'].is_zero is False:
            U,Y,UY=(clean(values[k]/values['Z0'])for k in MARKS[1:]);normalized={'U':U,'Y':Y,'UY':UY,'connected_UY':clean(UY-U*Y)}
        return row(IDS[0],tests,values=values,normalized=normalized,unnormalized_forms=source,
                   normalization_after_common_sum=True,whole_P_inserted_only_at_outer_contraction=True)
    return execute(IDS[0],impl,d)

def reflected_rp(d):
    def impl(d):
        words=need(d,'basis_words');Q0=matrix(need(d,'actual_Q0'))
        if not words or not Q0.rows:raise MissingInput('nonempty actual reflected word family')
        # This is an ORIGINAL word caller, not a supplied unrelated positive Gram.
        if callable(d.get('word_source')):
            theta=need(d,'theta_word');shift=need(d,'translate_word');T=need(d,'source_time')
            if not callable(theta)or not callable(shift):raise MissingInput('executable actual reflection/translation word maps')
            G=s.Matrix([[scalar(d['word_source'](theta(a),shift(b,T)))for b in words]for a in words])
            realization='actual_reflected_ordered_source_word_caller'
        elif d.get('ordered_word_amplitudes')is not None:
            amps=[matrix(x)for x in d['ordered_word_amplitudes']]
            norm=scalar(need(d,'reflection_normalization'))
            if norm.is_positive is not True:raise MissingInput('positive actual reflection normalization')
            if len(amps)!=len(words):raise ValueError('one actual sewn half-word amplitude per original word required')
            G=s.Matrix([[clean(s.trace(a.H*b)/norm)for b in amps]for a in amps])
            realization='actual_whole_sewn_halfword_amplitudes_no_semigroup_assumed'
        else:raise MissingInput('actual reflected word source or original sewn halfword amplitudes; input Gram alone is insufficient')
        tests={'word_count_matches_Q0':len(words)==Q0.rows,'actual_Q0_equals_reflected_word_form':equal(Q0,G),'actual_reflected_Q0_positive':psd(G)}
        qt,_=quotient(d,{'Q0':Q0});tests.update(qt)
        return row(IDS[1],tests,reflected_form=G,realization=realization,
                   no_Hamiltonian_or_reducing_block_hypothesis_added=True)
    return execute(IDS[1],impl,d)

def ward_sd(d):
    def impl(d):
        x=need(d,'coordinates');S=s.sympify(need(d,'action'));J=s.sympify(need(d,'coordinate_density'))
        v=[s.sympify(a)for a in need(d,'variation_vector')];O=s.sympify(need(d,'whole_marked_insertion'))
        integrate=need(d,'integrate_against_exp_minus_action');face=need(d,'integrate_boundary_face')
        bounds=need(d,'original_domain_bounds')
        if len(x)!=len(v)or len(x)!=len(bounds):raise ValueError('coordinate, variation and domain dimension mismatch')
        if not callable(integrate)or not callable(face):raise MissingInput('actual original bulk and boundary integration operations')
        # No log J or inverse fermion Pfaffian occurs; zeros of J are harmless.
        contact=sum(J*v[i]*s.diff(O,x[i])for i in range(len(x)))
        jacobian=sum(O*s.diff(J*v[i],x[i])for i in range(len(x)))
        action=-sum(O*J*v[i]*s.diff(S,x[i])for i in range(len(x)))
        inserted=clean(contact+jacobian+action)
        bulk=integrate(inserted)
        if bulk is None:raise MissingInput('actual original integrated Ward/SD insertion')
        flux=0
        for i,(lo,hi) in enumerate(bounds):
            for sign,endpoint in [(-1,lo),(1,hi)]:
                integrand=O*J*v[i]*s.exp(-S)
                actual=face(i,endpoint,integrand)
                if actual is None:raise MissingInput('actual original boundary/domain flux')
                flux+=sign*actual
        return row(IDS[2],{'integrated_variation_equals_original_boundary_flux':equal(scalar(bulk),scalar(flux))},
                   insertion_contact=contact,measure_jacobian_contact=jacobian,action_score_contact=action,
                   integrated_bulk=bulk,boundary_flux=flux,pointwise_EOM_filter_used=False)
    return execute(IDS[2],impl,d)

def whole_placement(d):
    def impl(d):
        B=matrix(need(d,'retained_embedding'))
        P=matrix(need(d,'whole_P'));cap=matrix(need(d,'whole_cap'))
        scalar_factor=scalar(need(d,'outer_scalar'));reference=scalar(need(d,'reference_scalar'))
        if reference.is_zero is not False:raise MissingInput('actual nonzero absolute reference scalar')
        if callable(d.get('whole_word_source')):
            recipes=need(d,'original_whole_word_recipes')
            full={k:matrix(d['whole_word_source'](need(recipes,k)))for k in MARKS}
            realization='actual_whole_word_caller_with_all_memory'
        else:
            K0,K1,K2=(matrix(a)for a in need(d,'original_interval_kernels'))
            U=matrix(need(d,'whole_U'));Y=matrix(need(d,'whole_Y'))
            full={'Z0':K2*K1*K0,'ZU':K2*K1*U*K0,'ZY':K2*Y*K1*K0,'ZUY':K2*Y*K1*U*K0}
            realization='actual_whole_interval_sewing_no_reducing_compression'
        actual=need(d,'observed_full_ordered_forms')
        tests={'whole_P_adjoint':equal(P,P.H),'whole_P_idempotent':equal(P*P,P)}
        projected={};values={}
        for name,A in full.items():
            tests[name+'_whole_sewing_matches']=equal(matrix(need(actual,name)),A)
            if d.get('observed_retained_source_forms')is not None:
                tests[name+'_actual_retained_packet_matches']=equal(matrix(need(d['observed_retained_source_forms'],name)),B.H*A*B)
            projected[name]=clean(B.H*P*A*P*B)
            values[name]=clean(scalar_factor*s.trace(cap*P*A*P)/reference)
        # Only the relations actually required by the ORIGINAL domain map are
        # tested. No universal commutation/reducing-block requirement is added.
        relations=need(d,'original_domain_transport_relations')
        for i,relation in enumerate(relations):tests['original_domain_transport_'+str(i)]=equal(matrix(need(relation,'actual')),matrix(need(relation,'required')))
        factors=need(d,'original_factor_transport_relations')
        for i,relation in enumerate(factors):tests['reference_eta_scalar_transport_'+str(i)]=equal(scalar(need(relation,'actual')),scalar(need(relation,'required')))
        return row(IDS[3],tests,whole_values=values,retained_projected_forms=projected,
                   realization=realization,
                   nonsemigroup_intervals_allowed=True,physical_projection_between_letters_used=False,
                   factor_relations_are_original_data_not_unit_defaults=True)
    return execute(IDS[3],impl,d)

def regulator_limit(d):
    def impl(d):
        kind=need(d,'kind')
        if kind=='exact_original_weighted_pushforward':
            fine=need(d,'fine_atoms');coarse=need(d,'coarse_atoms');p=need(d,'original_retained_map')
            if not fine or not coarse:raise MissingInput('nonempty actual fine and coarse original domains')
            if len({need(a,'id')for a in fine})!=len(fine)or len({need(a,'id')for a in coarse})!=len(coarse):raise ValueError('original atom identities must be unique')
            if not callable(p):raise MissingInput('actual original retained map')
            pushed={}
            for a in fine:
                k=p(a);w=scalar(need(a,'unnormalized_mass'));pushed[k]=pushed.get(k,0)+w
            target={need(a,'id'):scalar(need(a,'unnormalized_mass'))for a in coarse}
            keys=set(pushed)|set(target);tests={'full_unnormalized_pushforward_'+str(k):equal(pushed.get(k,0),target.get(k,0))for k in keys}
            result=row(IDS[4],tests,all_retained_bounded_words_covered=True,
                       scope_of_limit='this original finite weighted blocking map only; not ultraviolet tightness or all outer volumes')
            if d.get('required_limit_scope')!='finite_weighted_blocking':
                result.update(status='MISSING_INPUT',missing_inputs=['original matched infinite-regulator tail/tightness theorem and its numerical inputs'])
            return result
        if kind=='original_finite_transfer_cauchy_tail':
            T=matrix(need(d,'actual_transfer'));Pi=matrix(need(d,'actual_stationary_projector'));index=scalar(need(d,'n'))
            if not T.rows:raise MissingInput('nonempty actual original transfer')
            if index.is_integer is not True or index.is_nonnegative is not True:raise ValueError('exact nonnegative original iteration index required')
            n=int(index);r=matrix(T-Pi)
            eigen=list(r.eigenvals());q=max([abs(v)for v in eigen],default=s.Integer(0))
            tests={'transfer_adjoint':equal(T,T.H),'stationary_projector_adjoint':equal(Pi,Pi.H),'stationary_projector_idempotent':equal(Pi*Pi,Pi),
                   'stationary_left':equal(Pi*T,Pi),'stationary_right':equal(T*Pi,Pi),'strict_contraction':s.simplify(1-q).is_positive}
            if tests['strict_contraction']is not True:return row(IDS[4],tests)
            # Exact spectral theorem: ||R^k||<=q^k, so sum_(k>=n)R^k
            # has absolute operator norm <=q^n/(1-q), all subsequent indices.
            tail=clean(q**n/(1-q));exact=clean((s.eye(T.rows)-r).inv()*r**n)
            result=row(IDS[4],tests,q=q,n=n,tail_bound=tail,exact_remaining_series=exact,
                       all_remaining_iterations_analytically_bounded=True,
                       scope_of_limit='actual finite transfer iteration series, not a spatial ultraviolet field limit')
            if d.get('required_limit_scope')!='finite_transfer_iterations':
                result.update(status='MISSING_INPUT',missing_inputs=['original matched ultraviolet/volume limit bound beyond this finite transfer iteration theorem'])
            return result
        raise MissingInput('implemented original limit branch with actual weighted maps or analytic transfer tail; empirical N versus 2N is not a proof')
    return execute(IDS[4],impl,d)

FUNCTIONS=dict(zip(IDS,(common_four,reflected_rp,ward_sd,whole_placement,regulator_limit)))
def bind_current_caller(payload,result,entries,domain):
    """Bind independent mathematical inputs to the CURRENT actual atom result."""
    raw=need(domain,'whole_source_bundle')
    if isinstance(raw,(str,Path)):
        p=Path(raw);p=p if p.is_absolute()else Path(need(domain,'domain_directory'))/p
        bundle=json.loads(p.read_text())
    else:bundle=raw
    digest=hashlib.sha256(json.dumps(bundle,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
    if bundle.get('domain_id')!=need(domain,'domain_id')or bundle.get('native_inputs_sha256')!=need(domain,'_sha256'):
        raise ValueError('check bundle does not bind the actual runtime domain inputs')
    if need(payload,'source_bundle_sha256')!=digest:raise ValueError('check realization does not identify the actual current source bundle')
    if not entries:raise MissingInput('current original unnormalized contribution entries')
    atoms=[];expected={r['original_configuration_id']:r for r in need(bundle,'fiber_records')}
    if len(expected)!=len(bundle['fiber_records']):raise ValueError('original source atom identities must be unique')
    actual=[need(r,'configuration_id')for r in entries]
    if len(set(actual))!=len(actual)or set(actual)!=set(expected):raise ValueError('check requires the actual whole original atom coverage')
    for r in entries:
        original=expected[r['configuration_id']]
        if r.get('source_bundle_sha256')!=digest or r.get('domain_id')!=bundle['domain_id']or r.get('native_inputs_sha256')!=bundle['native_inputs_sha256']or r.get('original_measure_mass')!=original['original_measure_mass']:
            raise ValueError('actual contribution/source/domain binding mismatch')
        forms={z:need(need(need(r,'contributions'),w),'value')for z,w in zip(MARKS,('W0','WU','WY','WUY'))}
        atoms.append({'original_configuration_id':r['configuration_id'],'original_measure_mass':r['original_measure_mass'],
                      'unnormalized_scalar_weight':[1,0],'ordered_source_forms':forms})
    bound=dict(payload);bound['atoms']=atoms
    for name in ('gram_form','physical_selector_form','boundary_pairing_form','reference_scalar','scalar_factor'):
        bound[name]=need(bundle,name)
    bound['basis_words']=need(bundle,'retained_words')
    bound['observed_values']={k:need(result,k)for k in MARKS}
    summed={k:sum((scalar(a['original_measure_mass'])*matrix(a['ordered_source_forms'][k])for a in atoms),s.zeros(len(bound['basis_words'])))for k in MARKS}
    bound['actual_Q0']=summed['Z0'];bound['observed_retained_source_forms']=summed
    return bound
def ensemble_test_registry(domain):
    def callback(cid):
        def run(result,entries,actual_domain):
            payload=actual_domain.get('original_check_inputs',{}).get(cid,{})
            def checked(_):
                bound=bind_current_caller(payload,result,entries,actual_domain)
                checked_result=FUNCTIONS[cid](bound)
                checked_result['caller_binding_kind']=('actual_entries_and_observed_four_Z'if cid==IDS[0]else
                    'actual_aggregated_Q0_against_reflected_source_realization'if cid==IDS[1]else
                    'actual_retained_packets_against_whole_sewing'if cid==IDS[3]else
                    'source_input_digest_bound_supplied_original_integration_realization')
                checked_result['source_identity_not_proved_by_digest']=True
                if cid in (IDS[2],IDS[4]):
                    checked_result['applicability']='The supplied action/density/integrator or limit family must be the ORIGINAL same source. The implemented identity checks that realization; its hash alone does not prove that identification.'
                return checked_result
            return jsonable(execute(cid,checked,{}))
        return run
    return {cid:callback(cid)for cid in IDS}
def run_ensemble_checks(result,entries,domain):
    return [fn(result,entries,domain)for fn in ensemble_test_registry(domain).values()]
