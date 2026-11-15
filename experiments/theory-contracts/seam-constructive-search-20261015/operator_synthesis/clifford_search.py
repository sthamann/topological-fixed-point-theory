#!/usr/bin/env python3
"""Typed Majorana-word CEGIS recovery control, using exact Gaussian integers.

This recovers a known map in the GIVEN five-mode CAR source. It derives no
physical source, statistics, carrier interface, or P1/P2 selection principle.
"""
from __future__ import annotations
import argparse
from functools import lru_cache
from hashlib import sha256
import itertools
import json
from pathlib import Path
import random
import numpy as np

MODES=5
DIM=1<<MODES
HERE=Path(__file__).resolve().parent


def canonical(value):
    return json.dumps(value,sort_keys=True,separators=(',',':'))


def digest(value):
    return sha256(canonical(value).encode()).hexdigest()


def canonical_word(word):
    """gamma_i^2=1 and gamma_i gamma_j=-gamma_j gamma_i, including sign."""
    sign,mask=1,0
    for generator in word:
        if not isinstance(generator,int) or not 0<=generator<2*MODES:
            raise ValueError('Majorana generator outside the original ten generators')
        if (mask>>(generator+1)).bit_count()%2: sign=-sign
        mask^=1<<generator
    return sign,mask


def candidate_key(candidate):
    return str(candidate['support'])+':'+str(int(candidate['antilinear']))


def generate_candidates(seed=0):
    """Generate every word class by growing subsets, with neither parity seeded."""
    generators=list(range(2*MODES))
    random.Random(seed).shuffle(generators)
    for length in range(2*MODES+1):
        for word in itertools.combinations(generators,length):
            sign,support=canonical_word(word)
            for anti in (False,True):
                yield {'support':support,'antilinear':anti}


def pair(real,imag=None):
    real=np.asarray(real,dtype=np.int64)
    return real,np.zeros_like(real) if imag is None else np.asarray(imag,dtype=np.int64)


def adjoint(matrix):
    return matrix[0].T,-matrix[1].T


def multiply(a,b):
    return a[0]@b[0]-a[1]@b[1],a[0]@b[1]+a[1]@b[0]


def subtract(a,b):
    return a[0]-b[0],a[1]-b[1]


def phase(matrix,exponent):
    real,imag=matrix
    return ((real,imag),(-imag,real),(-real,-imag),(imag,-real))[exponent%4]


@lru_cache(None)
def annihilation(mode):
    real=np.zeros((DIM,DIM),dtype=np.int64)
    for state in range(DIM):
        if state&(1<<mode):
            real[state^(1<<mode),state]=(-1)**((state&((1<<mode)-1)).bit_count())
    return pair(real)


@lru_cache(None)
def word_action(support):
    if not isinstance(support,int) or not 0<=support<1<<(2*MODES):
        raise ValueError('Invalid ten-bit support')
    permutation=np.arange(DIM,dtype=np.int64)
    phases=np.zeros(DIM,dtype=np.int64)
    for generator in reversed(range(2*MODES)):
        if not support&(1<<generator): continue
        mode=generator//2
        for index,state in enumerate(permutation.copy()):
            state=int(state)
            exponent=2*((state&((1<<mode)-1)).bit_count()%2)
            if generator%2: exponent+=3 if state&(1<<mode) else 1
            phases[index]=(phases[index]+exponent)%4
            permutation[index]=state^(1<<mode)
    return permutation,phases


def transform(candidate,matrix):
    """U M U* or U conjugate(M) U*, entirely in integer real/imag blocks."""
    permutation,phases=word_action(candidate['support'])
    inverse=np.argsort(permutation)
    real=matrix[0][np.ix_(inverse,inverse)]
    imag=matrix[1][np.ix_(inverse,inverse)]
    if candidate['antilinear']: imag=-imag
    exponent=(phases[inverse,None]-phases[None,inverse])%4
    out_real=np.empty_like(real);out_imag=np.empty_like(imag)
    for k in range(4):
        selection=exponent==k
        r,i=phase((real,imag),k)
        out_real[selection]=r[selection];out_imag[selection]=i[selection]
    return out_real,out_imag


def fixture_key(fixture):
    return canonical(fixture)


def one_body_basis():
    for j in range(MODES): yield {'kind':'one_body','family':'diagonal','j':j,'k':j}
    for family in ('real_symmetric','imaginary_antisymmetric'):
        for j,k in itertools.combinations(range(MODES),2):
            yield {'kind':'one_body','family':family,'j':j,'k':k}


def complete_fixtures():
    return ([{'kind':'field','mode':j} for j in range(MODES)]
            +[{'kind':'charge'},{'kind':'quarter'}]+list(one_body_basis()))


def initial_fixtures():
    return ([{'kind':'field','mode':j} for j in range(MODES)]
            +[{'kind':'charge'},{'kind':'quarter'},
              {'kind':'one_body','family':'diagonal','j':0,'k':0},
              {'kind':'one_body','family':'real_symmetric','j':0,'k':1}])


@lru_cache(None)
def fixture_matrices(encoded,anti):
    fixture=json.loads(encoded)
    kind=fixture['kind']
    identity=pair(np.eye(DIM,dtype=np.int64))
    charge=pair(np.diag([state.bit_count() for state in range(DIM)]))
    if kind=='field':
        c=annihilation(fixture['mode'])
        return c,adjoint(c)
    if kind=='charge': return charge,subtract(pair(MODES*identity[0]),charge)
    if kind=='quarter':
        real=np.zeros((DIM,DIM),dtype=np.int64);imag=real.copy()
        for state in range(DIM):
            r,i=((1,0),(0,1),(-1,0),(0,-1))[state.bit_count()%4]
            real[state,state]=r;imag[state,state]=i
        quarter=(real,imag)
        # Antiunitarity conjugates i as well as the number operator.
        return quarter,phase(quarter,3) if anti else phase(adjoint(quarter),1)
    if kind!='one_body': raise ValueError('Unknown fixture kind')
    j,k=fixture['j'],fixture['k']
    first=multiply(adjoint(annihilation(j)),annihilation(k))
    trace=0
    if fixture['family']=='diagonal':
        if j!=k: raise ValueError('Diagonal basis fixture must have j=k')
        matrix=first;trace=1
    elif fixture['family']=='real_symmetric':
        other=adjoint(first)
        matrix=(first[0]+other[0],first[1]+other[1])
    elif fixture['family']=='imaginary_antisymmetric':
        matrix=phase(subtract(first,adjoint(first)),1)
    else: raise ValueError('Unknown Hermitian basis family')
    return matrix,subtract(pair(trace*identity[0]),matrix)


def counterexample(candidate,fixture):
    matrix,expected=fixture_matrices(fixture_key(fixture),candidate['antilinear'])
    actual=transform(candidate,matrix)
    for name,a,e in zip(('real','imaginary'),actual,expected):
        mismatch=np.argwhere(a!=e)
        if len(mismatch):
            row,column=map(int,mismatch[0])
            return {'fixture':fixture,'block':name,'row':row,'column':column,
                    'actual':int(a[row,column]),'expected':int(e[row,column])}
    return None


def verify_complete(candidate):
    for fixture in complete_fixtures():
        failure=counterexample(candidate,fixture)
        if failure: return failure
    return None


class ConstraintCache:
    def __init__(self,path=None):
        self.path=Path(path) if path else None
        self.state=json.loads(self.path.read_text()) if self.path and self.path.exists() else {}
        self.hits=0;self.invalidated=0

    def bind(self,fixtures):
        provenance={'fixtures_sha256':digest(fixtures),
                    'engine_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'modes':MODES}
        if self.state.get('provenance')!=provenance:
            self.invalidated+=len(self.state.get('entries',{}))
            self.state={'provenance':provenance,'entries':{}}
        return provenance['fixtures_sha256']

    def check(self,candidate,fixtures):
        self.bind(fixtures)
        key=candidate_key(candidate)
        if key in self.state['entries']:
            self.hits+=1
            return self.state['entries'][key]
        failure=next((failure for fixture in fixtures if (failure:=counterexample(candidate,fixture))),None)
        self.state['entries'][key]=failure
        return failure

    def save(self):
        if self.path: self.path.write_text(json.dumps(self.state,indent=2,sort_keys=True)+'\n')


def expression(candidate):
    word=[j for j in range(2*MODES) if candidate['support']&(1<<j)]
    return {'word':word,'support':candidate['support'],'antilinear':candidate['antilinear'],
            'expression':(' '.join('gamma_'+str(j) for j in word) or 'I')+(' K' if candidate['antilinear'] else ''),
            'global_sign':'canonical +1; conjugation-equivalent signs deduplicated'}


def search(run_dir,seed=0):
    run_dir=Path(run_dir);run_dir.mkdir(parents=True,exist_ok=True)
    state_file=run_dir/'constraints.json'
    previous=json.loads(state_file.read_text()) if state_file.exists() else {}
    allowed={fixture_key(f) for f in complete_fixtures()}
    fixtures=previous.get('fixtures',initial_fixtures())
    if any(fixture_key(f) not in allowed for f in fixtures): raise ValueError('Unknown persisted constraint')
    trace=list(previous.get('counterexample_trace',[]))
    reused=len(fixtures) if previous else 0
    candidates=list(generate_candidates(seed))
    cache=ConstraintCache(run_dir/'cache.json')
    rounds=[]
    for round_index in range(len(complete_fixtures())+1):
        fixture_hash=cache.bind(fixtures)
        survivors=[candidate for candidate in candidates if cache.check(candidate,fixtures) is None]
        rounds.append({'round':round_index,'fixture_hash':fixture_hash,
                       'active_constraints':len(fixtures),'survivors':[expression(c) for c in survivors]})
        learned=None
        for candidate in survivors:
            failure=verify_complete(candidate)
            if failure:
                learned=failure['fixture']
                trace.append({'candidate':expression(candidate),'counterexample':failure,
                              'old_fixture_hash':fixture_hash,'new_fixture_hash':digest(fixtures+[learned])})
                break
        if learned is None: break
        if learned in fixtures: raise ArithmeticError('Verifier returned an already enforced counterexample')
        fixtures.append(learned)
        cache.save()
        state_file.write_text(json.dumps({'fixtures':fixtures,'counterexample_trace':trace},indent=2,sort_keys=True)+'\n')
    else: raise ArithmeticError('CEGIS exceeded the complete finite basis')
    cache.save()
    state_file.write_text(json.dumps({'fixtures':fixtures,'counterexample_trace':trace},indent=2,sort_keys=True)+'\n')
    result={'outcome':'KNOWN_SOURCE_MAP_RECOVERY_CONTROL','seed':seed,
        'engine_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
        'candidate_count':len(candidates),'generated_from':'all word supports of lengths0..10, both linear/antilinear; no target word or matrix seeded',
        'arithmetic':'exact int64 Gaussian blocks; only bounded monomial conjugations and five-mode CAR bilinears, no float certification',
        'complete_basis':{'field_relations':5,'charge_relations':1,'quarter_relations':1,
                          'hermitian_one_body_diagonals':5,'real_symmetric':10,'imaginary_antisymmetric':10},
        'basis_completeness_argument':'Every Hermitian5x5 h is a real linear combination of these25 basis elements; both linear and antilinear conjugation are real-linear on this space.',
        'complete_fixture_sha256':digest(complete_fixtures()),
        'majorana_definition':'gamma_(2j)=c_j+c_j*, gamma_(2j+1)=i(c_j*-c_j)',
        'quarter_rules':{'linear':'J W J^-1 = i W^-1','antilinear':'J W J^-1 = -i W','W':'exp(i*pi*N/2)'},
        'rounds':rounds,'counterexample_trace':trace,'active_fixtures':fixtures,'fixture_sha256':digest(fixtures),
        'persisted_constraints_reused':reused,'cache_hits':cache.hits,'stale_cache_entries_invalidated':cache.invalidated,
        'survivors':[expression(c) for c in survivors],
        'source_premises':['given five-mode charged CAR carrier','given occupation basis and conjugation K',
                           'given number operator and all Hermitian one-body dynamics','given quarter operator exp(i*pi*N/2)'],
        'source_selection_proved':False,'new_physics_claim':False}
    (run_dir/'results.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir',type=Path,default=HERE/'runs/calibration')
    parser.add_argument('--seed',type=int,default=0)
    args=parser.parse_args()
    result=search(args.run_dir,args.seed)
    print(json.dumps({key:result[key] for key in ('outcome','candidate_count','survivors','fixture_sha256','stale_cache_entries_invalidated')},sort_keys=True))


if __name__=='__main__': main()
