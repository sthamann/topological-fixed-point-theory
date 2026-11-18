"""Exact words in a declared Clifford *-algebra freely joined with C[T].

T=-1 is selfadjoint, with no imposed relation to any Majorana or to T^2.
The source CAR interface is an assumption here, not inferred P1 physics.
"""
from __future__ import annotations
from dataclasses import dataclass
from functools import lru_cache
from collections.abc import Mapping

T=-1
SCHEMA_ID='tfpt.clifford-free-selfadjoint-word.v1'


def _check_m(m):
    if type(m) is not int or not 1<=m<=5:
        raise ValueError('m must be an integer from 1 through 5')


def _tokens(word,m):
    _check_m(m)
    word=tuple(word)
    if any(type(token) is not int or not (token==T or 0<=token<2*m) for token in word):
        raise ValueError('Word tokens must be T=-1 or original Majorana indices 0..2m-1')
    return word


def reduce_word(word,*,m=5):
    """Return (real phase, canonical tuple), reducing each T-separated segment."""
    word=_tokens(word,m)
    sign,mask=1,0
    out=[]
    for token in word:
        if token==T:
            out.extend(j for j in range(2*m) if mask&(1<<j))
            out.append(T)
            mask=0
        else:
            if (mask>>(token+1)).bit_count()%2:
                sign=-sign
            mask^=1<<token
    out.extend(j for j in range(2*m) if mask&(1<<j))
    return sign,tuple(out)


def star(word,*,m=5):
    """All generators are selfadjoint: reverse, then reduce, retaining phase."""
    return reduce_word(reversed(_tokens(word,m)),m=m)


def multiply(u,v,*,m=5):
    """Multiply unsigned words; callers multiply any existing term phases."""
    return reduce_word(_tokens(u,m)+_tokens(v,m),m=m)


def quarter(word,*,m=5):
    """Formal *-automorphism gamma_even->gamma_odd, gamma_odd->-gamma_even."""
    sign=1
    mapped=[]
    for token in _tokens(word,m):
        if token==T: mapped.append(T)
        elif token%2:
            mapped.append(token-1)
            sign=-sign
        else: mapped.append(token+1)
    reduced_sign,canonical=reduce_word(mapped,m=m)
    return sign*reduced_sign,canonical


@lru_cache(None)
def words(max_degree,*,m=5):
    """All distinct canonical basis words of total degree <= max_degree.

    Phases are coefficients, not duplicate basis elements. Sorted first by degree,
    then lexicographically; the empty tuple is the identity.
    """
    _check_m(m)
    if type(max_degree) is not int or max_degree<0:
        raise ValueError('max_degree must be a nonnegative integer')
    alphabet=(T,)+tuple(range(2*m))
    level={()}
    result=[()]
    for degree in range(1,max_degree+1):
        following=set()
        for prefix in level:
            for token in alphabet:
                _,word=reduce_word(prefix+(token,),m=m)
                if len(word)==degree: following.add(word)
        result.extend(sorted(following))
        level=following
    return tuple(result)


def word_to_json(word,*,m=5):
    """Sign-preserving JSON term; never discard a reduction phase."""
    phase,tokens=reduce_word(word,m=m)
    return {'phase':phase,'tokens':list(tokens)}


def word_from_json(value,*,m=5):
    """Read a signed JSON term, normalizing its tokens exactly."""
    if not isinstance(value,Mapping) or set(value)!={'phase','tokens'}:
        raise ValueError('A JSON term requires exactly phase and tokens')
    phase=value['phase']
    if type(phase) is not int or phase not in (-1,1) or not isinstance(value['tokens'],list):
        raise ValueError('JSON phase must be +1/-1 and tokens must be a list')
    sign,tokens=reduce_word(value['tokens'],m=m)
    return phase*sign,tokens


def schema(m=5):
    _check_m(m)
    return {'schema_id':SCHEMA_ID,'m':m,'T_token':T,
            'term':{'type':'object','required':['phase','tokens'],'additionalProperties':False,
                    'properties':{'phase':{'type':'integer','enum':[-1,1]},
                                  'tokens':{'type':'array','items':{'type':'integer','minimum':-1,'maximum':2*m-1}}}},
            'relations':['gamma_i*=gamma_i','gamma_i^2=1','gamma_i gamma_j=-gamma_j gamma_i for i!=j','T*=T'],
            'T_relations':'No commutator or T^2 reduction is imposed.',
            'quarter':'gamma_2j->gamma_(2j+1), gamma_(2j+1)->-gamma_2j, T->T',
            'source_scope':'Declared CAR interface and formal selfadjoint T, not a P1-derived physical source.'}


@dataclass(frozen=True)
class WordAlgebra:
    m:int=5

    def __post_init__(self): _check_m(self.m)
    def reduce_word(self,word): return reduce_word(word,m=self.m)
    def star(self,word): return star(word,m=self.m)
    def multiply(self,u,v): return multiply(u,v,m=self.m)
    def quarter(self,word): return quarter(word,m=self.m)
    def words(self,max_degree): return words(max_degree,m=self.m)
    def word_to_json(self,word): return word_to_json(word,m=self.m)
    def word_from_json(self,value): return word_from_json(value,m=self.m)
    def schema(self): return schema(self.m)
