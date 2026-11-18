"""Scoped all-word charge obstruction and original finite source controls.

The obstruction assumes the explicitly stated primitive field inventory and
symmetry action. It is not a universal no-go for the original P1 source.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
VERIFICATION = ROOT/'verification'


def need(condition, label, checks):
    if not bool(condition):
        raise ArithmeticError('Failed: '+label)
    checks.append(label)


def load_original(name):
    path = VERIFICATION/(name+'.py')
    sys.path.insert(0, str(VERIFICATION))
    spec = importlib.util.spec_from_file_location('_source_'+name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, {'path':str(path.relative_to(ROOT)),
                    'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}


def add(a,b):
    return ((a[0]+b[0])%2,(a[1]+b[1])%2)


def assess_expression(expr):
    """Assess a JSON expression under the declared source inventory.

    Nodes use {op: primitive, name: ...}, {op: proposed_sector_insertion,
    source_ref: ...}, {op: product|sum, args: [...]}, or
    {op: adjoint|derivative, arg: ...}. This only checks tags and types.
    """
    primitive_tags = {'matter_spinor': (1, 1), 'tensor_boson': (0, 0),
                      'scalar_quarter_current': (0, 0)}
    count = 0

    def visit(node, depth):
        nonlocal count
        count += 1
        if count > 100 or depth > 24 or type(node) is not dict:
            raise ValueError('invalid expression or complexity limit')
        op = node.get('op')
        if type(op) is not str:
            raise ValueError('operation must be a string')
        if op == 'primitive':
            if set(node) != {'op', 'name'} or type(node['name']) is not str:
                raise ValueError('malformed primitive')
            name = node['name']
            if name == 'clifford_vector':
                raise ValueError('Clifford vector is an internal operator, not a local field')
            if name not in primitive_tags:
                raise ValueError('unknown primitive')
            return {primitive_tags[name]}, True
        if op == 'proposed_sector_insertion':
            if (set(node) != {'op', 'source_ref'} or type(node['source_ref']) is not str
                    or not node['source_ref'].strip()):
                raise ValueError('proposed sector requires a nonempty source_ref string')
            # Its type is not established by its name. Overapproximate all
            # four tags so composition cannot create a false exclusion.
            return {(0, 0), (0, 1), (1, 0), (1, 1)}, False
        if op in ('adjoint', 'derivative'):
            if set(node) != {'op', 'arg'}:
                raise ValueError('malformed unary expression')
            return visit(node['arg'], depth + 1)
        if op in ('product', 'sum'):
            if set(node) != {'op', 'args'} or type(node['args']) is not list or not node['args']:
                raise ValueError('malformed composite expression')
            tags, bound = ({(0, 0)} if op == 'product' else set()), True
            for item in node['args']:
                child_tags, child_bound = visit(item, depth + 1)
                bound = bound and child_bound
                tags = ({add(a, b) for a in tags for b in child_tags}
                        if op == 'product' else tags | child_tags)
            return tags, bound
        raise ValueError('unknown operation')

    try:
        tags, bound = visit(expr, 1)
    except (ValueError, RecursionError) as error:
        return {'possible_target_tag': False, 'source_binding_verified': False,
                'verdict': 'typing_error', 'reason': str(error),
                'candidate_physical_proof': False}
    possible = (0, 1) in tags
    return {'possible_target_tag': possible, 'source_binding_verified': False,
            'tags': [list(tag) for tag in sorted(tags)],
            'sector_insertion_status': None if bound else 'unverified_original_binding',
            'verdict': ('excluded_scoped' if not possible else 'needs_original_source'),
            'candidate_physical_proof': False}


def pauli_clifford_10():
    """Exact 0/±1/±i Clifford action on the 32-dimensional spinor space."""
    identity=np.eye(2,dtype=complex)
    x=np.array([[0,1],[1,0]],complex)
    y=np.array([[0,-1j],[1j,0]],complex)
    z=np.diag([1,-1]).astype(complex)
    output=[]
    for j in range(5):
        for a in (x,y):
            factors=[z if k<j else a if k==j else identity for k in range(5)]
            matrix=factors[0]
            for factor in factors[1:]: matrix=np.kron(matrix,factor)
            output.append(matrix)
    return output


def check():
    checks=[]
    # A word's pair of symmetry eigenvalues is the additive homomorphism
    # from the free *-word monoid to Z2 x Z2; conjugation leaves each tag.
    generators=[(1,1),(0,0)]
    reachable={(0,0)}
    while True:
        enlarged=reachable | {add(x,g) for x in reachable for g in generators}
        if enlarged==reachable: break
        reachable=enlarged
    need(reachable=={(0,0),(1,1)},'primitive parity and central tag closure',checks)
    need((0,1) not in reachable,'bosonic spinor current tag excluded in this inventory',checks)
    need(add((1,1),(1,1))==(0,0),'two spinor source fields form neutral tag',checks)
    enlarged=reachable | {add(x,(0,1)) for x in reachable}
    need((0,1) in enlarged and len(enlarged)==4,'explicit disorder primitive removes this tag obstruction',checks)
    need(all(add(x,(0,0))==x for x in reachable),'neutral source insertions preserve tags',checks)
    need(all(add(x,x)==(0,0) for x in reachable),'adjoint has the same Z2 tags',checks)

    # Finite exact matrix realization: a 2π spin rotation is -I on spinors;
    # conjugation by it is +I on the Clifford vector.  This is not gamma_F.
    gammas=pauli_clifford_10(); eye32=np.eye(32,dtype=complex)
    need(all(np.array_equal(g@g,eye32) for g in gammas),'ten Clifford generators square to identity',checks)
    need(all(np.array_equal(gammas[i]@gammas[j]+gammas[j]@gammas[i],np.zeros((32,32)))
             for i in range(10) for j in range(i)), 'all ten Clifford generators anticommute',checks)
    z2=(gammas[0]@gammas[1])@(gammas[0]@gammas[1])
    need(np.array_equal(z2,-eye32),'central spin rotation acts as -I on spinors',checks)
    need(all(np.array_equal(z2@g@z2.conj().T,g) for g in gammas),
         'same central element conjugates Clifford vector trivially',checks)

    original, source = load_original('v252_full_finite_triple')
    # Integer complex entries remain exact under these small matrix products.
    y=np.array([[1,2j,0],[0,2,1],[1,0,3]],complex)
    mr=np.array([[1,0,1j],[0,2,0],[1j,0,3]],complex)
    D=original.full_D(y,y.conj(),y,y.conj(),mr)
    J=original.Kswap(); gamma=original.grading(); central=-np.eye(original.NH)
    need(original.NH==96 and D.shape==(96,96),'original v252 full 96 carrier',checks)
    need(np.array_equal(central@D,D@central),'central spin action commutes with full v252 D',checks)
    need(np.array_equal(central@J,J@central.conjugate()),'central spin action commutes with antilinear J',checks)
    need(not np.array_equal(central,gamma),'central spin action is not finite chirality gamma_F',checks)

    # Actual original E8 root generator; D8 integer and half-spinor classes.
    roots, roots_source=load_original('v1_e8_glue')
    values=roots.e8_roots()
    half=[r for r in values if np.all(np.abs(r)==0.5)]
    ints=[r for r in values if np.all(r==np.rint(r))]
    spin16_4=[r for r in half if sum(v<0 for v in r[:5])%2==0]
    need(len(ints)==112 and len(half)==128 and len(values)==240
         and all(int(2*r[0])%2==0 for r in ints)
         and all(int(2*r[0])%2==1 for r in half),
         'original E8 roots split 112 plus 128 with correct central-square signs',checks)
    need(len(spin16_4)==64 and len(half)-len(spin16_4)==64,
         'marked D5 spinor and conjugate each carry 64 roots',checks)
    need(len(ints)+8==120,'D8 even roots plus Cartan give 120',checks)

    # Ward illustration: invariant vacuum; charged field one-point vanishes,
    # neutral positive two-point can be nonzero.
    rotation=np.diag([1,1j]); vacuum=np.array([1,0],complex)
    X=np.array([[0,0],[1,0]],complex)
    one=vacuum.conj()@X@vacuum
    two=vacuum.conj()@X.conj().T@X@vacuum
    need(np.array_equal(rotation@vacuum,vacuum) and np.array_equal(rotation@X@rotation.conj().T,1j*X),
         'explicit invariant state and charge-one operator',checks)
    need(one==0 and two==1,'charged one-point zero but neutral two-point positive',checks)

    return {'schema':'tfpt.source-charge-obstruction.v1','status':'pass','checks':checks,
            'source_paths':[source,roots_source],
            'algebraic_argument':'The commuting parity and central involutions have product Pz². Every listed primitive is fixed by Pz², so all sums, products and adjoints of them are fixed; the fixed algebra is norm closed for a continuous symmetry action. A target of tag (0,1) changes sign and cannot be a nonzero element of that algebra. Equivariant derivations and neutral insertions preserve the conclusion when defined on that algebra.',
            'tag_class':[[0,0],[1,1]],'target_tag':[0,1],
            'explicit_escape_tag':[0,1],
            'ward_scope':'Invariant state and nontrivial character force a charged one-point to zero; neutral correlators need not vanish. The finite matrix example illustrates this identity and is not the physical TFPT state.',
            'source_scope':'Given primitive Grassmann Spin10 spinors and tensorial bosons, plus equivariant operations. A disorder/sector primitive may change this inventory. E8 current is bosonic; 16 internal spinor states are not 16 Majorana generators. The Spin action z²=-I96 is supplied on the v252 carrier, and commutation with D and J is checked; v252 alone does not derive that representation. z² on H_F is distinct from gamma_F and field parity.',
            'not_proved':['P1 source has only this inventory','all possible disorder sectors are absent',
                          'physical 3+1D matter, E8 current source or global P1 no-go'],
            'physical_source_identified':False,'fundamental_solution':False}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path)
    claims = {
        'matter-word-spinor-current': 'even_matter_words_contain_bosonic_spinor_current',
        'charged-one-point': 'invariant_state_has_nonzero_charge_one_onepoint',
        'clifford-vector-spinor-typing': 'clifford_vector_is_spinor_current_under_same_spin_action',
    }
    parser.add_argument('--check', choices=sorted(claims))
    args=parser.parse_args()
    result=check()
    if args.check:
        result['id'] = args.check
        result['inference_claims'] = [{'id': claims[args.check], 'holds': False}]
    encoded=json.dumps(result,indent=2)+'\n'
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        with args.out.open('x') as stream: stream.write(encoded)
    print(encoded,end='')


if __name__=='__main__': main()
