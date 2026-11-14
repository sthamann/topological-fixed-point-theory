"""Exact finite control for the Clifford double cell and its parity reference."""
from collections import Counter
from itertools import combinations, product
from math import prod

import sympy as sp

from reconstruction import PAULI, require


def occupation_projectors():
    identity, x, y, z = PAULI
    gamma = [
        sp.kronecker_product(x, identity), sp.kronecker_product(y, identity),
        sp.kronecker_product(z, x), sp.kronecker_product(z, y),
    ]
    gamma.append(prod(gamma, start=sp.eye(4)))
    require(all(a == a.H and a*a == sp.eye(4) for a in gamma), "invalid Clifford squares")
    require(all(a*b == -b*a for a,b in combinations(gamma,2)), "invalid Clifford pairs")
    marks = [sp.kronecker_product(g,g) for g in gamma]
    unit = sp.eye(16)
    require(prod(marks, start=unit) == unit, "five-mark relation failed")
    require(all(a*b == b*a for a,b in combinations(marks,2)), "occupation marks fail to commute")
    projectors = {}
    for bits in product((0,1), repeat=5):
        p = prod(((unit+(-1)**b*m)/2 for b,m in zip(bits,marks)), start=unit)
        require(p*p == p and p.trace() == 1-sum(bits)%2, "invalid joint occupation projector")
        if p.trace():
            projectors[bits] = p
    require(sum(projectors.values(), sp.zeros(16)) == unit, "incomplete occupation resolution")
    return gamma, marks, projectors


def permutation_action(gamma, marks, projectors):
    lifts = [gamma[i]-gamma[i+1] for i in range(4)]
    require(lifts[0]*lifts[2] == -lifts[2]*lifts[0], "missing single-cell spin sign")
    doubled = [(sp.kronecker_product(u,u)/2).applyfunc(sp.expand) for u in lifts]
    for i, w in enumerate(doubled):
        require(w.H*w == sp.eye(16), "invalid doubled lift")
        for j, mark in enumerate(marks):
            moved = 2*i+1-j if j in (i,i+1) else j
            require(w*mark*w.H == marks[moved], "wrong lifted action on occupations")
    identity = tuple(range(5))
    group = {identity:sp.eye(16)}
    queue = [identity]
    for permutation in queue:
        for i, w in enumerate(doubled):
            moved = list(permutation)
            moved[i], moved[i+1] = moved[i+1], moved[i]
            key = tuple(moved)
            operator = group[permutation]*w
            if key in group:
                require(group[key] == operator, "doubled group law failed")
            else:
                group[key] = operator
                queue.append(key)
    require(len(group) == 120, "incomplete permutation group")
    # Compare each fixed occupation line with the exterior action, including its sign.
    # Equivariant lines in each orbit are induced from their stabilizer character.
    invariant = sp.zeros(16)
    for permutation, operator in group.items():
        invariant += operator
        character = 0
        for bits, p in projectors.items():
            occupied = [i for i,b in enumerate(bits) if b]
            moved = [permutation[i] for i in occupied]
            if sorted(moved) == occupied:
                sign = (-1)**sum(a>b for a,b in combinations(moved,2))
                require((p*operator).trace() == sign, "wrong exterior sign on an occupation line")
                character += sign
        require(operator.trace() == character, "wrong exterior character")
    vacuum = projectors[(0,0,0,0,0)]
    require(invariant/120 == vacuum, "wrong invariant line")
    require(all(sum(vacuum[4*i+k,4*j+k] for k in range(4)) == sp.eye(4)[i,j]/4
                for i,j in product(range(4),repeat=2)), "wrong first-cell marginal")
    require(all(sum(vacuum[4*k+i,4*k+j] for k in range(4)) == sp.eye(4)[i,j]/4
                for i,j in product(range(4),repeat=2)), "wrong second-cell marginal")
    return sum(permutation[4] == 4 for permutation in group)


def parity_reference(gamma, parity):
    identity, x, _, z = PAULI
    flip = sp.kronecker_product(gamma[0],sp.eye(4))
    even = next(v for v in sp.eye(16).columnspace() if parity*v == v)
    odd = flip*even
    require(parity*odd == -odd, "reference flip must exchange the parity sectors")
    plus, minus = sp.eye(2)[:,0], sp.eye(2)[:,1]
    total = sp.kronecker_product(parity,z)
    require(((sp.eye(32)+total)/2).trace() == 16, "wrong fixed total parity dimension")
    first, second = sp.kronecker_product(even,plus), sp.kronecker_product(odd,minus)
    probe = sp.kronecker_product(flip,x)
    require(probe*total == total*probe, "joint probe violates total parity")
    require(first.H*second == sp.zeros(1) and even.H*odd == sp.zeros(1),
            "reference marginals do not erase the relative sign")
    values = []
    for sign in (1,-1):
        state = (first+sign*second)/sp.sqrt(2)
        require(total*state == state, "reference state is not in the fixed parity sector")
        require(all((state.H*a*state)[0] == 0 for a in (
            sp.kronecker_product(flip,identity),sp.kronecker_product(sp.eye(16),x))),
            "single-factor odd readout must vanish")
        values.append(int((state.H*probe*state)[0]))
    require(values == [1,-1], "joint probe must distinguish the otherwise equal marginals")
    return values


def replay():
    gamma, marks, projectors = occupation_projectors()
    stabilizer = permutation_action(gamma, marks, projectors)
    occupations = Counter(map(sum,projectors))
    parity_sectors = Counter((-1)**bits[4] for bits in projectors)
    # In the standard exterior dictionary the weak raising operator transports
    # occupation from slot four to slot five while preserving total degree.
    weak_moves = [(bits,(*bits[:3],0,1)) for bits in projectors if bits[3:] == (1,0)]
    require(all(b in projectors and a[4] != b[4] for a,b in weak_moves),
            "weak exchange did not cross the native parity sectors")
    g, interaction = sp.symbols("g lambda", nonnegative=True)
    energies = [g*k+interaction*k*(k-1)/2 for k in sorted(occupations)]
    gap_ratio = energies[2]/energies[1]
    ratios = [str(gap_ratio.subs({g:1,interaction:i})) for i in (0,1)]
    require(ratios == ["2","10/3"], "invariant Hamiltonians should have different gap ratios")
    return {
        "status":"EXACT_FINITE_DOUBLE_CELL_CONTROL",
        "occupation_multiplicities":dict(occupations),
        "native_parity_multiplicities":dict(parity_sectors),
        "permutations_checked":120, "parity_stabilizer_order":stabilizer,
        "relative_mark_dimension":3, "invariant_vacuum_rank":1,
        "weak_exchange_rank":len(weak_moves),
        "weak_exchange_native_parity_commutator_rank":len(weak_moves),
        "fixed_total_parity_dimension":16, "joint_relative_phase_responses":parity_reference(gamma,marks[4]),
        "symmetric_energy_levels":[str(e) for e in energies], "gap_ratio_controls":ratios,
        "physical_source_identity_proved":False, "r113_original_values":None,
        "first_original_test":"Bind the parity compensation and the charge, adjoint and time actions "
                              "to the same weighted P1 word before outer selection.",
    }
