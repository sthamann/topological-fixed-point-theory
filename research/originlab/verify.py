"""Independent exact controls and discriminators for OriginLab."""
from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import gcd, lcm
from pathlib import Path
import gzip
import hashlib
import json
import sympy as sp

from language import positive

CAR = Path("experiments/theory-contracts/seam-positive-kernel-search-20261015/runs/source-clustering")
CAR_ARCHIVE = Path(__file__).with_name("fixtures")/"r006"
CAR_ORIGINAL_HASHES = {
    "moments.json":"7429db68e991b7c5351361f607ac9fb2dfb62e3c99b6945813fdb5160f5bda9c",
    "min-trace.json":"b713dd585daeb2984a25c9c4bd4bee4f5cef0a06895ea7b54e5f14f5ad0949ed",
    "max-trace.json":"9809d92e023a3011ef6d418766d2f6b10e37d35430177b1e6f2cdcdc2c81a577",
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_text())


def r102(low, high):
    # The presentation supplies the D and A dual norm residues, not an answer.
    solutions = []
    visited = 0
    for n in range(low, high+1):
        norms = (0, 1, Fraction(n+1, 4), Fraction(n+1, 4))
        for b in range(1, n):
            if gcd(n, b) != 1:
                continue
            for index, norm in enumerate(norms):
                visited += 1
                total = norm+Fraction((n-1)*b*b, n)
                if total.denominator == 1:
                    d_order = 1 if index == 0 else 2 if index == 1 or (n+1)%2 == 0 else 4
                    order = lcm(d_order, n)
                    solutions.append(dict(n=n, charge=b, dual_class=index, norm=str(total),
                                          index=order, determinant=str(Fraction(4*n, order**2))))
    return dict(status="exact_bounded", bounds=[low, high], checked=visited, solutions=solutions,
                conclusion="Within the standard pairing, integrality implies n | 4b²; primitive b implies n | 4.",
                premises=["standard positive D/A dual norm residues", "integral generator", "primitive A class"])


def read_matrix(data):
    def number(x):
        return sp.Rational(*x["re"]) + sp.I*sp.Rational(*x["im"])
    return sp.SparseMatrix([[number(x) for x in row] for row in data])


def scalar_recovery(value):
    rational = Fraction(value).limit_denominator(64)
    if abs(float(rational)-value) > 1e-6:
        raise ValueError("stored decimal has no admitted rational recovery")
    return sp.Rational(rational.numerator, rational.denominator)


def car_candidate(data, gap, floor):
    record = data["exact_matrices"]
    state, transfer, quarter = (read_matrix(record[k]) for k in ("Omega", "T", "quarter"))
    letters = [read_matrix(g) for g in record["gammas"]]
    identity = sp.eye(transfer.rows)
    zero = sp.zeros(transfer.rows)
    if state.H*state != sp.ones(1) or transfer*state != state or quarter*state != state:
        raise ValueError("invalid normalized stationary state")
    if quarter.H*quarter != identity or quarter**4 != identity or quarter*transfer != transfer*quarter:
        raise ValueError("invalid compatible quarter action")
    for i, a in enumerate(letters):
        paired = letters[i ^ 1] * (1 if i % 2 == 0 else -1)
        if a.H != a or a*a != identity or quarter*a*quarter.H != paired:
            raise ValueError("invalid Clifford generator or quarter action")
        if any(a*b+b*a != zero for b in letters[:i]):
            raise ValueError("invalid CAR anticommutation")
    if not positive(transfer-floor*identity) or not positive((1-gap)*identity+gap*state*state.H-transfer):
        raise ValueError("transfer fails exact positivity or declared clustering bound")
    return state, {-1: transfer, **dict(enumerate(letters))}


def separator(left, right, max_length, max_questions=20):
    """Exact reachable pair-space search. Linear dependence safely prunes suffixes."""
    tokens = sorted(left[1])
    if set(tokens) != set(right[1]):
        raise ValueError("candidate alphabets differ")
    left_size = left[0].rows
    row = left[0].H.row_join(-right[0].H)
    basis = {}
    scanned = Counter()
    found = []
    layer = [((), left[0].col_join(right[0]))]
    for degree in range(max_length+1):
        following = []
        for word, vector in layer:
            scanned[degree] += 1
            if (row*vector)[0] != 0:
                first = (left[0].H*vector[:left_size, :])[0]
                second = (right[0].H*vector[left_size:, :])[0]
                found.append(dict(word=list(word), responses=[str(first), str(second)]))
                if len(found) == max_questions:
                    break
            reduced = list(vector)
            for index, pivot in sorted(basis.items()):
                coefficient = reduced[index]
                if coefficient:
                    reduced = [sp.expand(a-coefficient*b) for a,b in zip(reduced,pivot)]
            index = next((i for i,x in enumerate(reduced) if x != 0), None)
            if index is None:
                continue
            coefficient = reduced[index]
            basis[index] = [sp.expand(x/coefficient) for x in reduced]
            if degree < max_length:
                for token in tokens:
                    v = left[1][token]*vector[:left_size, :]
                    w = right[1][token]*vector[left_size:, :]
                    following.append(((token, *word), v.col_join(w)))
        if found or not following:
            break
        layer = following
    return dict(questions=found, shortest_length=len(found[0]["word"]) if found else None,
                nodes_examined=dict(scanned), reachable_rank=len(basis),
                equivalence_proved=not found and not following and degree < max_length)


def car_replay(root, max_length):
    folder = root/CAR
    records,hashes,archives = {},{},[]
    for name,expected in CAR_ORIGINAL_HASHES.items():
        path = folder/name
        if path.is_file():
            raw = path.read_bytes()
        else:
            raw = gzip.decompress((CAR_ARCHIVE/(name+".gz")).read_bytes())
            if hashlib.sha256(raw).hexdigest()!=expected:
                raise ValueError(f"altered archived R006 fixture {name}")
            archives.append(name)
        records[name] = json.loads(raw)
        hashes[name] = hashlib.sha256(raw).hexdigest()
    source = records["moments.json"]
    gap, floor = (scalar_recovery(source[k]) for k in ("cluster_gap", "transfer_floor"))
    candidates = [car_candidate(records[name],gap,floor) for name in ("min-trace.json","max-trace.json")]
    for state, ops in candidates:
        cache = {(): state}

        def vector(word):
            if word not in cache:
                cache[word] = ops[word[0]]*vector(word[1:])
            return cache[word]

        for row in source["moments"]:
            expected = scalar_recovery(row["re"])+sp.I*scalar_recovery(row["im"])
            if (state.H*vector(tuple(row["word"])))[0] != expected:
                raise ValueError(f"moment mismatch for {row['word']}")
    return dict(status="exact_finite_controls", shared_moments=len(source["moments"]),
                evidence_inputs=hashes,archived_inputs=archives,
                discriminator=separator(*candidates, max_length),
                premises=["declared finite CAR representation", "rational recovery of supplied moments"])


def clock_test():
    z = sp.Symbol("z")
    # Permutation on four coordinates, then quotient the invariant constant line.
    cycle = sp.zeros(4)
    for index in range(4):
        cycle[(index+1)%4, index] = 1
    cartan = sp.div(cycle.charpoly(z).as_expr(), z-1, z)[0]
    spatial = sp.Matrix([[0,-1,0],[1,0,0],[0,0,1]]).charpoly(z).as_expr()
    return dict(status="refuted" if cartan != spatial else "undecided",
                cartan_polynomial=str(sp.factor(cartan)), spatial_polynomial=str(sp.factor(spatial)),
                scope="direct similarity of the two specified quarter actions")


def eigenvalue_test():
    # The ten balanced 3+3 lines in the sum-zero subspace of C^6.
    # The resulting projectors are closed under coordinate permutations.
    labels = [subset for subset in combinations(range(6),3) if 0 in subset]
    states = [sp.Matrix([1 if j in label else -1 for j in range(6)]) for label in labels]
    projectors = [v*v.T/6 for v in states]
    carrier = sp.eye(6)-sp.ones(6)/6
    if sum(projectors,sp.zeros(6)) != 2*carrier:
        raise ArithmeticError("Bell projectors do not resolve the carrier")
    gram = sp.Matrix([[(a*b).trace()/2 for b in projectors] for a in projectors])
    if gram.eigenvals() != {sp.Integer(1):1,sp.Rational(4,9):9}:
        raise ArithmeticError("unexpected nonzero spectrum of the Bell return")
    # Convex mixture of two CP trace preserving maps on the carrier:
    # K(X)=sum Pi X Pi/2 and D(X)=Tr(X) I_carrier/5.
    # D vanishes on the traceless subspace. K and D share covariance and vacuum.
    weights = [sp.Integer(0),sp.Rational(1,4),sp.Rational(1,2)]
    return dict(status="nonunique", mixture_weights=list(map(str,weights)),
                witnesses=[str((1-w)*sp.Rational(4,9)) for w in weights],
                bell_projectors=len(projectors), bell_rank=5,
                scope="CPTP S6-covariant positive selfadjoint maps on the same Bell carrier, same maximally mixed stationary state and superoperator rank; no requirement that every map be the unchanged original instrument",
                exact_construction="K_w=(1-w)K+wD, D(X)=Tr(X)I/5")
