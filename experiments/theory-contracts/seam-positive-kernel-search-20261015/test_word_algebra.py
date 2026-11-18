"""Independent exact Fock checks, using arbitrary noncommuting Hermitian T."""
import importlib.util
import itertools
import json
from pathlib import Path
import random
import unittest
import sys
import sympy as s

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('word_algebra',HERE/'word_algebra.py')
a=importlib.util.module_from_spec(spec);sys.modules[spec.name]=a;spec.loader.exec_module(a)


def fock(m,seed):
    identity=s.eye(2)
    x=s.Matrix([[0,1],[1,0]])
    y=s.Matrix([[0,-s.I],[s.I,0]])
    z=s.diag(1,-1)
    matrices={}
    for j in range(m):
        for parity,pauli in enumerate((x,y)):
            factors=[z]*j+[pauli]+[identity]*(m-j-1)
            matrices[2*j+parity]=s.kronecker_product(*factors)
    rng=random.Random(seed)
    dim=2**m
    t=s.zeros(dim)
    for j in range(dim):
        t[j,j]=rng.randint(-2,2)
        for k in range(j):
            value=rng.randint(-2,2)+s.I*rng.randint(-2,2)
            t[j,k]=value;t[k,j]=s.conjugate(value)
    matrices[a.T]=t
    return matrices


def evaluate(word,matrices):
    result=s.eye(next(iter(matrices.values())).rows)
    for token in word: result=result*matrices[token]
    return s.expand(result)


def term_eval(term,matrices):
    phase,word=term
    return phase*evaluate(word,matrices)


class WordAlgebraTests(unittest.TestCase):
    def test_CAR_segments_and_T_are_distinct(self):
        self.assertEqual(a.reduce_word([2,1,2,1]),(-1,()))
        self.assertEqual(a.reduce_word([2,1,-1,2,1]),(1,(1,2,-1,1,2)))
        self.assertEqual(a.reduce_word([0,-1,0]),(1,(0,-1,0)))
        self.assertEqual(a.reduce_word([-1,-1]),(1,(-1,-1)))
        self.assertNotEqual(a.reduce_word([0,-1]),a.reduce_word([-1,0]))
        self.assertNotEqual(a.reduce_word([0,-1,0]),a.reduce_word([-1]))
        matrices=fock(2,123)
        self.assertNotEqual(evaluate([0,-1],matrices),evaluate([-1,0],matrices))

    def test_exhaustive_small_words_against_independent_Fock(self):
        matrices=fock(2,29)
        for length in range(5):
            for word in itertools.product((-1,0,1,2,3),repeat=length):
                self.assertEqual(term_eval(a.reduce_word(word,m=2),matrices),evaluate(word,matrices))

    def test_associativity_star_and_quarter_exactly(self):
        rng=random.Random(11)
        for m in range(1,6):
            alphabet=(-1,)+tuple(range(2*m))
            for _ in range(100):
                u,v,w=[tuple(rng.choice(alphabet) for _ in range(rng.randrange(7))) for _ in range(3)]
                p,uv=a.multiply(u,v,m=m);q,left=a.multiply(uv,w,m=m)
                r,vw=a.multiply(v,w,m=m);t,right=a.multiply(u,vw,m=m)
                self.assertEqual((p*q,left),(r*t,right))
                p,us=a.star(u,m=m);q,uss=a.star(us,m=m)
                self.assertEqual((p*q,uss),a.reduce_word(u,m=m))
                p,uv=a.multiply(u,v,m=m);q,lhs=a.star(uv,m=m)
                r,vs=a.star(v,m=m);t,us=a.star(u,m=m);h,rhs=a.multiply(vs,us,m=m)
                self.assertEqual((p*q,lhs),(r*t*h,rhs))
                phase,word=1,u
                for _ in range(4):
                    sign,word=a.quarter(word,m=m);phase*=sign
                self.assertEqual((phase,word),a.reduce_word(u,m=m))
                p,uv=a.multiply(u,v,m=m);q,lhs=a.quarter(uv,m=m)
                r,uq=a.quarter(u,m=m);t,vq=a.quarter(v,m=m);h,rhs=a.multiply(uq,vq,m=m)
                self.assertEqual((p*q,lhs),(r*t*h,rhs))

    def test_star_and_quarter_against_independent_matrix_substitution(self):
        rng=random.Random(41)
        for m in (1,2,3):
            matrices=fock(m,m+40)
            mapped={-1:matrices[-1]}
            for j in range(2*m): mapped[j]=matrices[j+1] if j%2==0 else -matrices[j-1]
            for _ in range(30):
                word=tuple(rng.choice(tuple(matrices)) for _ in range(rng.randrange(7)))
                self.assertEqual(term_eval(a.star(word,m=m),matrices),evaluate(word,matrices).adjoint())
                self.assertEqual(term_eval(a.quarter(word,m=m),matrices),evaluate(word,mapped))

    def test_enumeration_is_complete_unique_and_not_phase_duplicated(self):
        for m in range(1,4):
            for degree in range(4):
                brute={a.reduce_word(word,m=m)[1] for length in range(degree+1)
                       for word in itertools.product((-1,)+tuple(range(2*m)),repeat=length)}
                enumerated=a.words(degree,m=m)
                self.assertEqual(set(enumerated),brute)
                self.assertEqual(len(enumerated),len(set(enumerated)))
                self.assertTrue(all(a.reduce_word(word,m=m)==(1,word) for word in enumerated))
        self.assertEqual(len(a.words(2,m=5)),78)
        self.assertEqual(len(a.words(3,m=5)),419)
        self.assertEqual(len(a.words(4,m=5)),2245)
        products={a.multiply(u,v,m=5)[1] for u in a.words(2,m=5) for v in a.words(2,m=5)}
        self.assertEqual(products,set(a.words(4,m=5)))

    def test_JSON_schema_and_configured_API(self):
        algebra=a.WordAlgebra(2)
        payload=algebra.word_to_json([2,1,-1,3,3])
        self.assertEqual(payload,{'phase':-1,'tokens':[1,2,-1]})
        self.assertEqual(algebra.word_from_json(json.loads(json.dumps(payload))),(-1,(1,2,-1)))
        self.assertEqual(algebra.schema()['m'],2)
        self.assertEqual(algebra.quarter([-1]),(1,(-1,)))
        for bad in (0,6,True):
            with self.assertRaises(ValueError): a.WordAlgebra(bad)
        for word in ([-2],[4],[True]):
            with self.assertRaises(ValueError): algebra.reduce_word(word)
        with self.assertRaises(ValueError): algebra.word_from_json({'phase':0,'tokens':[]})
        with self.assertRaises(ValueError): algebra.words(-1)


if __name__=='__main__': unittest.main()
