"""Exact whole-family comparison against implementation-bound source providers.

A candidate may propose an intertwiner and variable names. It cannot register
an original provider or turn a hash/role/status into physical authority. The
fixed-frame transport is the existing L20261102G same-source theorem: transport
all operators, reference factors, domains and time together. No semigroup is
assumed for an arbitrary original interval family.
"""
from dataclasses import dataclass, field, replace
from typing import Callable
import sympy as s

IDS = ('primitive_original_source_identity', 'original_marked_time_word_match')
ORIGINAL_SCOPE = 'primitive_original_source'
PROVIDER_LOCATION = (
    'experiments/theory-contracts/original-paper-source-reconciliation-20261102/'
    'contribution_20261106Q/runtime/original_source_adapter.py:'
    '_evaluate_contribution_unchecked -> evaluate_contribution; '
    'upstream whole source/time family, not its four supplied jet matrices')


class SourceFamilyMissing(ValueError):
    def __init__(self, message, details=None):
        super().__init__(message)
        self.details = details or {}


@dataclass(frozen=True)
class SourceTimeFamily:
    """Entire represented kernel, with code-owned complete source coordinates.

    ``carried_operators`` includes the actual generator, U/Y, reference/domain
    operators or reflection whenever they are part of this original caller.
    No new universal operator is required by this comparison. Source and time
    symbol maps describe the complete parameters of this represented family;
    they are not a selected list of four Taylor coefficients.
    Typed extra legs act on these same coordinates: antilinear coefficients
    use W J conjugate(W)^-1, bilinear forms use congruence, and domain relations
    are pulled back. A one-particle form is never put on a different Fock
    space merely because a dimension matches. ``reference_factors`` contains
    the actually paired factors including their measure/orientation, not an
    isolated Pfaffian whose compensating Jacobian was discarded.
    """
    kernel: object
    source_symbols: dict
    time_symbols: dict
    P: object
    cap: object
    reference_scalar: object
    outer_scalar: object
    carried_operators: dict = field(default_factory=dict)
    reference_factors: dict = field(default_factory=dict)
    antilinear_operators: dict = field(default_factory=dict)
    bilinear_forms: dict = field(default_factory=dict)
    domain_relations: dict = field(default_factory=dict)


@dataclass(frozen=True)
class SourceProvider:
    """Host-registered producer. Never constructed from candidate JSON."""
    load: Callable
    acceptance_scope: str
    location: str


def equal(a, b):
    a, b = (s.Matrix(a), s.Matrix(b)) if isinstance(a, s.MatrixBase) else (s.sympify(a), s.sympify(b))
    if isinstance(a, s.MatrixBase) and a.shape != b.shape:
        return False
    difference = a-b
    entries = list(difference) if isinstance(difference, s.MatrixBase) else [difference]
    def identity(x):
        expression=s.simplify(x)
        answer=expression.is_zero
        if answer is not None:
            return answer
        # A nonzero polynomial is not an identical family, even though it
        # vanishes at selected parameter values. This is exact coefficient
        # comparison, not floating point sampling or a random identity test.
        if expression.free_symbols:
            try:
                coefficients=s.Poly(expression,*sorted(expression.free_symbols,key=str)).coeffs()
                signs=[s.simplify(v).is_zero for v in coefficients]
                if False in signs:
                    return False
                if all(v is True for v in signs):
                    return True
            except s.PolynomialError:
                pass
        return None
    tests = [identity(x) for x in entries]
    return False if False in tests else True if all(x is True for x in tests) else None


def verdict(tests):
    return 'FAIL' if False in tests.values() else 'NOT_PROVED' if None in tests.values() else 'PASS'


def _missing(reason, location=PROVIDER_LOCATION, **extra):
    return {key: {'status': 'MISSING_INPUT', 'missing_inputs': [reason],
                  'provider_location': location, **extra} for key in IDS}


def _renaming(candidate_map, original_map, kind):
    if not original_map or len(set(original_map.values()))!=len(original_map):
        raise SourceFamilyMissing('nonempty independently parameterized original '+kind+' family')
    if not isinstance(candidate_map, dict) or set(candidate_map) != set(original_map):
        raise SourceFamilyMissing('complete original '+kind+' symbol correspondence')
    result = {}
    for name, target in original_map.items():
        proposed = candidate_map[name]
        if not isinstance(proposed, s.Symbol) or not isinstance(target, s.Symbol):
            raise ValueError(kind+' correspondence allows symbol renaming, not a new scale or nonlinear reparameterization')
        if proposed in result:
            raise ValueError('distinct original '+kind+' coordinates were merged')
        result[proposed] = target
    return result


def transport_family(original, unitary):
    """Existing full fixed-frame U theorem, retaining same time and scalars.

    This transforms a supplied family; it does not produce the primitive
    source. For a field-dependent/moving frame the original connection and
    derivative terms must instead be supplied by its owning producer.
    """
    W = s.Matrix(unitary);K = s.Matrix(original.kernel)
    if W.shape != K.shape or equal(W.H*W, s.eye(K.rows)) is not True or equal(W*W.H, s.eye(K.rows)) is not True:
        raise ValueError('the whole fixed frame must be exactly unitary')
    parameters = set(original.source_symbols.values()) | set(original.time_symbols.values())
    if any(x.has(*parameters) for x in W) if parameters else False:
        raise ValueError('this existing transport is fixed; moving-frame source/time contacts cannot be discarded')
    move = lambda A: W*s.Matrix(A)*W.H
    return replace(original, kernel=move(K), P=move(original.P), cap=move(original.cap),
                   carried_operators={k: move(v) for k, v in original.carried_operators.items()},
                   antilinear_operators={k: W*s.Matrix(v)*W.T for k,v in original.antilinear_operators.items()},
                   bilinear_forms={k: W.conjugate()*s.Matrix(v)*W.H for k,v in original.bilinear_forms.items()},
                   domain_relations={k: s.Matrix(v)*W.H for k,v in original.domain_relations.items()})


def compare_whole_family(candidate, family, packet):
    """Compare exact functions, all their source orders, and all original legs."""
    if not isinstance(family, SourceTimeFamily):
        raise SourceFamilyMissing('implementation-produced SourceTimeFamily')
    if 'matrix' not in candidate or 'dictionary' not in candidate:
        raise SourceFamilyMissing('constructed whole kernel and outer dictionary')
    if not isinstance(packet, dict):
        raise SourceFamilyMissing('proposed whole source/time intertwiner, not a four-jet match')
    sources = _renaming(packet.get('source_symbols'), family.source_symbols, 'source')
    times = _renaming(packet.get('time_symbols'), family.time_symbols, 'time')
    if set(sources) & set(times):
        raise ValueError('a source variable was aliased to physical time')
    if 'intertwiner' not in packet:
        raise SourceFamilyMissing('whole Hilbert-space intertwiner after the original null quotient')
    W=s.Matrix(packet['intertwiner']);K=s.Matrix(candidate['matrix']);O=s.Matrix(family.kernel)
    if K.rows != K.cols or O.rows != O.cols or W.shape != (K.rows, O.rows):
        raise ValueError('whole source/intertwiner dimensions do not match')
    mapping={**sources, **times};parameters=set(family.source_symbols.values()) | set(family.time_symbols.values())
    # A proposed frame depending on the varied sources/time needs the actual
    # original moving-frame contacts, not this fixed-frame identity.
    # A renamed candidate coordinate must not hide a frame depending on its
    # original target coordinate. Both belong to this varied source/time
    # family; the existing theorem here covers fixed frames only.
    varied_parameters = set(mapping) | parameters
    if any(x.has(*varied_parameters) for x in W) if varied_parameters else False:
        raise SourceFamilyMissing('actual moving-frame source/time contacts; the existing fixed-frame theorem alone is insufficient')
    rename=lambda value:s.Matrix(value).subs(mapping, simultaneous=True)
    dic=candidate['dictionary']
    needed={'P','cap','reference_scalar','outer_scalar'}
    if needed-set(dic):
        raise SourceFamilyMissing('complete original P/cap/reference/absolute dictionary')
    source_tests={
        'whole_isometry':equal(W.H*W,s.eye(O.rows)),
        'whole_onto_not_a_nonreducing_compression':equal(W*W.H,s.eye(K.rows)),
        'entire_source_kernel_all_orders':equal(rename(K)*W,W*O),
        'whole_P':equal(rename(dic['P'])*W,W*s.Matrix(family.P)),
        'whole_cap':equal(rename(dic['cap'])*W,W*s.Matrix(family.cap)),
        'absolute_reference':equal(s.sympify(dic['reference_scalar']).subs(mapping,simultaneous=True),family.reference_scalar),
        'outer_scalar':equal(s.sympify(dic['outer_scalar']).subs(mapping,simultaneous=True),family.outer_scalar),
    }
    proposed_operators=packet.get('carried_operators',{})
    if set(proposed_operators)!=set(family.carried_operators):
        raise SourceFamilyMissing('all actually present original time/mark/domain/reference operators')
    for name,A in family.carried_operators.items():
        source_tests['operator:'+name]=equal(rename(proposed_operators[name])*W,W*s.Matrix(A))
    for field_name,rule in (
        ('antilinear_operators',lambda A,B: equal(A*W.conjugate(),W*B)),
        ('bilinear_forms',lambda A,B: equal(W.T*A*W,B)),
        ('domain_relations',lambda A,B: equal(A*W,B)),
    ):
        expected=getattr(family,field_name);proposed=packet.get(field_name,{})
        if set(proposed)!=set(expected):
            raise SourceFamilyMissing('all actually present original '+field_name)
        for name,A in expected.items():
            source_tests[field_name+':'+name]=rule(rename(proposed[name]),s.Matrix(A))
    proposed_factors=packet.get('reference_factors',{})
    if set(proposed_factors)!=set(family.reference_factors):
        raise SourceFamilyMissing('all actually present paired reference/eta/orientation factors')
    for name,value in family.reference_factors.items():
        source_tests['factor:'+name]=equal(s.sympify(proposed_factors[name]).subs(mapping,simultaneous=True),value)
    # Exact equality as a function includes every time, not selected samples.
    time_tests={'whole_time_family_with_same_parameter':source_tests['entire_source_kernel_all_orders']}
    for name,A in family.carried_operators.items():
        time_tests['actual_carried_operator:'+name]=source_tests['operator:'+name]
    # All scalar/reference/operator legs are needed by marked words as well.
    time_tests.update({'same_source_transport:'+k:v for k,v in source_tests.items()})
    return {
        IDS[0]:{'status':verdict(source_tests),'tests':source_tests,
                'scope':'whole represented source family; equality of all source orders, not four jets'},
        IDS[1]:{'status':verdict(time_tests),'tests':time_tests,
                'scope':'same original source variables, whole time family and carried marks; no semigroup inferred'},
    }


def original_identity_results(candidate, providers=None, context=None, required_scope=ORIGINAL_SCOPE):
    """Only implementation-owned providers can earn original-source credit."""
    registry={} if providers is None else providers
    packet=candidate.get('original_packet')
    provider_id=packet.get('provider_id','original.common_source') if isinstance(packet,dict) else 'original.common_source'
    provider=registry.get(provider_id)
    if not isinstance(provider,SourceProvider) or not callable(provider.load):
        return _missing('code-bound original whole source/time provider; a source ID/hash or four supplied matrices is not that producer')
    if provider.acceptance_scope!=required_scope:
        return _missing('provider of the required original scope',provider.location,
                        supplied_provider_scope=provider.acceptance_scope,required_scope=required_scope)
    try:
        family=provider.load(context)
    except SourceFamilyMissing as exc:
        return _missing(str(exc),provider.location,**exc.details)
    except Exception as exc:
        # A host-bound producer failure is not a counterexample to the
        # candidate. Keep the producer and comparison failure axes separate.
        return {key:{'status':'ERROR','error_stage':'original_provider_load',
                     'error_type':type(exc).__name__,'error':str(exc),
                     'provider_location':provider.location}for key in IDS}
    try:
        results=compare_whole_family(candidate,family,packet)
    except SourceFamilyMissing as exc:
        return _missing(str(exc),provider.location,**exc.details)
    except (ValueError,TypeError,KeyError,s.ShapeError)as exc:
        return {key:{'status':'FAIL','counterexample':str(exc),'provider_location':provider.location}for key in IDS}
    for result in results.values():
        result.update(provider_location=provider.location,acceptance_scope=provider.acceptance_scope,
                      hash_or_role_used_as_proof=False,candidate_self_attestation_used=False)
    return results


def providers_for_repository(repo):
    """Use the actual installed original contribution path, never candidate code.

    The current installed body consumes four supplied forms. Calling it is
    useful evidence, but it cannot reconstruct a whole time/source family
    from those four jets. A subsequently derived full original producer must
    be registered by the implementation, not selected by a candidate label.
    """
    from pathlib import Path
    import hashlib
    import importlib.util
    import json
    base=Path(repo)/'experiments/theory-contracts/original-paper-source-reconciliation-20261102/contribution_20261106Q'
    adapter_path=base/'runtime/original_source_adapter.py'
    domain_path=base/'domain/native-diagnostic-v1.json'

    def load(context):
        context=context or {}
        binding=context.get('runtime_binding') or {}
        if not adapter_path.is_file():
            raise SourceFamilyMissing('installed original source contribution operation',
                                      {'missing_original_file':str(adapter_path.relative_to(Path(repo)))})
        spec=importlib.util.spec_from_file_location('_r_original_source_identity_native',adapter_path)
        native=importlib.util.module_from_spec(spec);spec.loader.exec_module(native)
        domain=binding.get('domain')
        if domain is None:
            if not domain_path.is_file():
                raise SourceFamilyMissing('actual original source domain for the installed contribution call')
            raw=domain_path.read_bytes();domain=json.loads(raw)
            domain.update(_sha256=hashlib.sha256(raw).hexdigest(),domain_directory=str(domain_path.parent))
        plan=native.integration_plan(domain)
        items=plan.get('items') or []
        probe=items[0]['configuration'] if items else {
            'schema':'tfpt.original-source-unit.v1','original_configuration_id':None}
        assembled={'candidate_id':probe.get('native_candidate_id'),
                   'domain_id':domain.get('domain_id'),'native_inputs_sha256':domain.get('_sha256'),
                   'h_N_role':'not_aliased_to_original_joint_generator'}
        attempt=native.evaluate_contribution(probe,assembled,domain)
        details={'actual_native_call':{
            'operation':'original_source_adapter.evaluate_contribution',
            'status':attempt.get('status'),
            'missing_operations':attempt.get('missing_operations',[]),
            'four_contributions_present':all(attempt.get('contributions',{}).get(k) is not None
                                             for k in ('W0','WU','WY','WUY'))},
            'existing_transport':'L20261102G fixed full frame preserves each supplied whole S/Z/OS source; it is not a new source producer'}
        raise SourceFamilyMissing(
            'actual complete original source/time family upstream of the installed contribution operation; '
            'four supplied jet forms do not determine all source orders and histories',details)

    return {'original.common_source':SourceProvider(load,ORIGINAL_SCOPE,PROVIDER_LOCATION)}
