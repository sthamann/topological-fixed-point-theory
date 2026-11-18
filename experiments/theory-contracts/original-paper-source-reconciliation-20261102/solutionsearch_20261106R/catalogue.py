"""Curated construction grammar. Every new relation is explicitly a hypothesis.

The complete even retained-kernel basis is available, not only the quartic
template. Coefficient programs use an explicit ring of original constants,
algebraic units and formal source variables. These relations do not acquire
physical provenance by being syntactically well typed.
"""
from pathlib import Path
import hashlib
import sympy as s
from search.synthesis import InputBinding, Unavailable, digest, evidence_json
from source_recovery import native_constructors as native


def fp(value):
    return digest(evidence_json(value))


def inputs():
    ju, jy = s.symbols('jU jY', real=True)
    result = {}
    def add(name, typ, value):
        result[name] = InputBinding(typ, value, fp(value))
    add('source_symbols', 'SourceSymbolPair', (ju, jy))
    add('original_endpoint_marks', 'SourceMarkPair', Unavailable(
        'Actual original U/Y endpoint operators with their common physical time binding have not been supplied to this retained kernel chart. A labelled hypothesis is not their derivation.'))
    add('original_source_dictionary', 'WholeSourceDictionary', Unavailable(
        'The actual original P/cap/reference/absolute/time/marks on this retained construction remain to be computed; no unit physical defaults.'))
    for mask in range(256):
        if mask.bit_count() % 2 == 0:
            add('even_mask_'+str(mask), 'EvenGrassmannMask', mask)
    for name, expr in [('zero', s.S.Zero), ('one', s.S.One), ('jU', ju), ('jY', jy)]:
        add(name, 'SourceScalar', {'expression': expr, 'hypotheses': [],
                                 'provenance': [{'role': 'formal algebra unit or source symbol', 'name': name}]})
    # These are explicit retained-chart SEARCH HYPOTHESES. They are never
    # labelled as the original physical mode/field dictionary.
    for n in (2, 4):
        chart = {'modes': n, 'mode_chart_id': 'hypothesis-retained-kernel-chart-'+str(n),
                 'hypotheses': ['Identification of this retained chart with the original whole marked source remains to be derived.']}
        add('chart_'+str(n), 'ModeChart', chart)
        roles = {name: {'role': 'hypothesis', 'description': text} for name, text in {
            'P': 'Candidate identity selector on the retained chart, not the derived original physical selector.',
            'cap': 'Candidate ordinary coherent trace cap, not an inferred primitive cap or state.',
            'reference': 'Candidate unit relative reference; original absolute reference must be matched.',
            'absolute': 'Candidate unit outer scalar; the searched body remains inside the kernel.',
            'time': 'Candidate whole marked kernel, with no inferred semigroup or physical time generator.',
            'U': 'The jU derivative is a proposed insertion, not automatically the original Higgs norm.',
            'Y': 'The jY derivative is a proposed insertion, not automatically the local HYPF Hamiltonmark.',
        }.items()}
        add('dictionary_'+str(n), 'WholeSourceDictionary',
            {'P': s.eye(2**n), 'cap': s.eye(2**n), 'reference_scalar': s.S.One,
             'outer_scalar': s.S.One, 'bindings': roles})
        # All even basis monomials. At this finite interface arbitrary even
        # kernels over the declared coefficient ring can be constructed.
        for mask in range(1 << (2*n)):
            if mask.bit_count() % 2 == 0:
                add('monomial_'+str(n)+'_'+str(mask), 'EvenMonomial',
                    {**chart, 'mask': mask, 'source_symbols': (ju, jy)})
    return result


def kernel_term(monomial, coefficient):
    return {**monomial, 'coefficients': {monomial['mask']: coefficient['expression']},
            'hypotheses': monomial.get('hypotheses', [])+
                coefficient.get('hypotheses', [])+
                ['This weighted kernel coefficient relation is a search hypothesis.'],
            'provenance': coefficient.get('provenance', [])}


def kernel_add(left, right):
    for key in ('modes', 'mode_chart_id', 'source_symbols'):
        if left[key] != right[key]:
            raise ValueError('Kernel addition requires the same labelled source chart')
    coefficients = dict(left['coefficients'])
    for mask, coefficient in right['coefficients'].items():
        coefficients[mask] = s.expand(coefficients.get(mask, 0)+coefficient)
    return {**left, 'coefficients': coefficients,
            'hypotheses': left.get('hypotheses', [])+right.get('hypotheses', []),
            'provenance': left.get('provenance', [])+right.get('provenance', [])}


def scalar_negate(value):
    return {**value, 'expression': -value['expression']}


def scalar_inverse(value):
    e = value['expression']
    if e.is_zero is not False:
        raise ValueError('Inverse needs an actually nonzero coefficient, not an assumed denominator')
    return {**value, 'expression': 1/e}


def action_polynomial(program):
    return {**program, 'coefficients': dict(program['nil_action'])}


def kernel_of_operator(operator):
    return {**operator, 'coefficients': native.native_module().kernel(operator['matrix'])}


def operations():
    v = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    extra = [
        ('construct.kernel_term', ('EvenMonomial', 'SourceScalar'), 'EvenGrassmannKernel', kernel_term),
        ('construct.kernel_add', ('EvenGrassmannKernel', 'EvenGrassmannKernel'), 'EvenGrassmannKernel', kernel_add),
        ('construct.negate', ('SourceScalar',), 'SourceScalar', scalar_negate),
        ('construct.inverse_nonzero', ('SourceScalar',), 'SourceScalar', scalar_inverse),
        ('construct.action_polynomial', ('KernelActionProgram',), 'EvenGrassmannKernel', action_polynomial),
        ('construct.kernel_of_operator', ('EvenSourceOperator',), 'EvenGrassmannKernel', kernel_of_operator),
    ]
    return native.registry_operations()+[
        {'id': i, 'input_types': it, 'output_type': ot, 'execute': fn, 'version': v}
        for i, it, ot, fn in extra]


def leaf(name):
    return {'op': 'input', 'id': name, 'args': []}


def apply(name, *args):
    return {'op': 'apply', 'id': name, 'args': list(args)}


def starting_programs():
    """Productive programs evaluated before enumerating shorter fragments."""
    for n in (2, 4):
        for body, q in [(apply('native.c3'), apply('native.c6_ratio')),
                        (leaf('one'), apply('construct.negate', apply('native.c6_ratio')))]:
            action = apply('native.quartic_source_program', leaf('chart_'+str(n)),
                           body, q, leaf('source_symbols'))
            yield apply('native.source_jets', apply('native.nil_action', action),
                        leaf('dictionary_'+str(n)))
