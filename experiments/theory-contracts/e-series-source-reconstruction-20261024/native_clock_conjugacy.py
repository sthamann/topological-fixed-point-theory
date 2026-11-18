#!/usr/bin/env python3
"""Complete E8 toral-four-torsion test for the ACTUAL native field clock.

This exact certificate plus the compact-group conjugacy theorem proves an
inner VOA conjugacy, not an originally selected microscopic conjugator.
"""
import argparse
from collections import Counter
import contextlib
import io
from itertools import product
import json
from pathlib import Path
import runpy

import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent


def need(value, label):
    if not bool(value):
        raise RuntimeError(label)


def bits(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length()-1
        mask ^= bit


def bitset(values):
    return sum(1 << i for i in values)


def calculate():
    with contextlib.redirect_stdout(io.StringIO()):
        old = runpy.run_path(str(HERE.parent/"e-series-joint-source-20261023/joint_clock.py"))
    B, Bi = old["B"], old["Bi"]
    gram = B.T*B/4
    need(gram.det() == 1, "unimodular E8 integral torus coordinates")
    need(all(gram[i,i] == 2 for i in range(8)) and
         all(gram[i,j] in (0,-1) for i in range(8) for j in range(8) if i != j),
         "native basis is an actual simply-laced simple-root basis")
    roots = [tuple(int(x) for x in v) for v in old["roots"]]
    root_lookup = {v:i for i,v in enumerate(roots)}
    root_coords = sp.Matrix(old["coordinates"].T.tolist())
    pairings = np.array((gram*root_coords).T.tolist(),dtype=np.int16)
    need(np.max(np.abs(pairings)) <= 2, "int16 product bound <=48, no overflow")
    candidates = np.array(list(product(range(4),repeat=8)),dtype=np.int16)
    phases = (candidates@pairings.T) % 4
    characters = np.stack([np.sum(phases == k,axis=1) for k in range(4)],axis=1)
    characters[:,0] += 8
    target_character = np.array([60,64,60,64])
    selected_indices = np.flatnonzero(np.all(characters == target_character,axis=1))
    selected = {tuple(map(int,row)) for row in candidates[selected_indices]}
    need(len(candidates) == 65536 and len(selected) == 15120, "complete toral census")

    v = sp.Matrix([sp.Rational(1,2)]*5+[0]*3)
    # a=4v in L; doubled standard coordinates are 8v.
    representative = tuple(int(x) % 4 for x in Bi*(8*v))
    need(representative in selected, "actual carrier parity lift kappa")
    need(representative == (0,1,3,2,0,2,0,0), "ordered original-basis representative")
    integer_gram = [[int(gram[j,k]) for k in range(8)] for j in range(8)]
    permutations = []
    for j in range(8):
        beta = tuple(int(x) for x in B[:,j])
        permutation = []
        for alpha in roots:
            d = sum(x*y for x,y in zip(alpha,beta))
            need(d % 4 == 0, "root reflection coefficient")
            permutation.append(root_lookup[tuple(x-(d//4)*y for x,y in zip(alpha,beta))])
        permutations.append(permutation)
    seed_a3 = bitset(i for i,r in enumerate(roots)
                    if sum(x != 0 for x in r) == 2 and all(x == 0 for x in r[:5]))
    need(seed_a3.bit_count() == 12, "original A3 normal frame")
    orbit = {representative:seed_a3}
    queue = [representative]
    for point in queue:
        for j in range(8):
            image = list(point)
            image[j] = (image[j]-sum(integer_gram[j][k]*point[k] for k in range(8))) % 4
            image = tuple(image)
            a3 = bitset(permutations[j][i] for i in bits(orbit[point]))
            if image not in orbit:
                orbit[image] = a3
                queue.append(image)
            else:
                need(orbit[image] == a3, "clock stabilizer preserves actual A3 factor")
    need(set(orbit) == selected, "ALL matching toral clocks form ONE actual Weyl orbit")
    root_kernels = {}
    for row_id in selected_indices:
        point = tuple(map(int,candidates[row_id]))
        mask = int.from_bytes(np.packbits(phases[row_id] == 0,bitorder="little").tobytes(),"little")
        root_kernels.setdefault(mask,[]).append(point)
    need(len(root_kernels) == 7560 and {len(a) for a in root_kernels.values()} == {2},
         "two primitive clocks per full compiler flag")
    dot_matrix = np.array(roots,dtype=np.int16)@np.array(roots,dtype=np.int16).T
    orthogonal = [bitset(np.flatnonzero(row == 0).tolist()) for row in dot_matrix]
    a3_count = Counter()
    for kernel,pair in root_kernels.items():
        need(tuple((-x) % 4 for x in pair[0]) == pair[1], "the two clocks are inverses")
        a3 = orbit[pair[0]]
        need(a3 == orbit[pair[1]], "inverse clocks determine the same A3 frame")
        core = (1 << 240)-1
        for i in bits(a3):
            core &= orthogonal[i]
        need(core.bit_count() == 40 and (a3 | core) == kernel,
             "actual fixed root kernel is precisely D5 core plus A3 normal")
        a3_count[a3] += 1
    need(len(a3_count) == 7560 and set(a3_count.values()) == {1},
         "bijection of unoriented primitive clocks and all tetrahedron frames")

    def degree(r):
        return sum(r[:5]) % 4

    kappa_counts = Counter(degree(r) for r in roots)
    kappa_counts[0] += 8
    need([kappa_counts[i] for i in range(4)] == [60,64,60,64], "full kappa field spectrum")
    for r in roots:
        need((2*degree(r)) % 4 == 2*(r[0] % 2), "kappa square is original E8/D8 spinor sign")
    need(old["power_basis_phases"](old["J"],0,4) == (1,)*8, "actual Jhat order four on full fields")

    # Retain the WHOLE native cyclic compiler, not just its ninth power.
    # Root cycle phases use the actual representative over c, not a fitted
    # torus spectrum. Its 4th/9th powers are the inherited sigmahat/Jhat.
    cstd = old["J_standard"]*old["sigma_standard"]
    c_character = [0]*12
    unseen = set(range(240))
    c_mask = min(row["c_input_character"] for row in old["solutions"])
    c_cycles = Counter()
    while unseen:
        start = min(unseen)
        point, cycle = start, []
        while point not in cycle:
            cycle.append(point)
            point = root_lookup[tuple(cstd*sp.Matrix(roots[point]))]
        need(point == start and len(cycle) in (4,12), "actual compiler root cycles")
        unseen -= set(cycle)
        c_cycles[len(cycle)] += 1
        phase, _ = old["power_phase_at"](old["c"],c_mask,len(cycle),
                                         tuple(int(x) for x in old["coordinates"][start]))
        need(phase == 1, "literal inherited compiler cycle phase")
        for k in range(0,12,12//len(cycle)):
            c_character[k] += 1
    x = sp.Symbol("x")
    need(sp.factor(old["c"].charpoly(x).as_expr()) == (x*x+1)**2*(x**4-x*x+1),
         "actual Cartan compiler character")
    for k in (3,9):
        c_character[k] += 2
    for k in (1,5,7,11):
        c_character[k] += 1
    need(c_cycles == {12:19,4:3} and c_character == [22,20,19,24,19,20,22,20,19,24,19,20],
         "complete actual C12 field character")
    sigma_candidates = np.array(list(product(range(3),repeat=8)),dtype=np.int16)
    kappa_degrees = np.array([degree(r) for r in roots],dtype=np.int16)
    compiler_phases = (3*kappa_degrees + 4*((sigma_candidates@pairings.T) % 3)) % 12
    compiler_characters = np.stack([np.sum(compiler_phases == k,axis=1) for k in range(12)],axis=1)
    compiler_characters[:,0] += 8
    matching = {tuple(map(int,a)) for a in sigma_candidates[
        np.all(compiler_characters == np.array(c_character),axis=1)]}
    need(len(sigma_candidates) == 6561 and len(matching) == 480, "complete cubic choice census")
    carrier_roots = [i for i,r in enumerate(roots) if degree(r) == 0]
    reflection_data = [(tuple(int(t) for t in root_coords[:,i]),
                        tuple(int(t) for t in pairings[i])) for i in carrier_roots]
    sigma_representative = min(matching)
    sigma_orbit, sigma_queue = {sigma_representative}, [sigma_representative]
    for a in sigma_queue:
        for vector,covector in reflection_data:
            d = sum(u*v for u,v in zip(a,covector))
            image = tuple((u-d*v) % 3 for u,v in zip(a,vector))
            if image not in sigma_orbit:
                sigma_orbit.add(image)
                sigma_queue.append(image)
    need(sigma_orbit == matching, "ALL matching compiler clocks form one actual carrier-Weyl orbit")
    sigma_torus = (B/2)*sp.Matrix(sigma_representative)/3
    # Keeping the old sigma frozen is expressly impossible: this crossing
    # root has opposite kappa phases before and after the old sigma.
    crossing = (0,0,0,0,2,2,0,0)
    need(degree(crossing) != degree(tuple(old["sigma_standard"]*sp.Matrix(crossing))),
         "a joint conjugator must transport sigma, not freeze it")
    result = {
        "research_id":"MATH.SEAM.ESERIES.SOURCE-RECONSTRUCTION.06.NATIVE-CLOCK",
        "finite_certificate":"PASS",
        "toral_four_torsion_candidates":65536,
        "points_with_actual_native_J_character":len(selected),
        "matching_actual_Weyl_orbits":1,
        "kappa_integral_torus_residue":representative,
        "kappa_half_coordinate_generator":["1/2"]*5+["0"]*3,
        "kappa_square":"gamma: +1 on D8, -1 on its spinor sector",
        "V1_character_1_i_minus1_minusi":[60,64,60,64],
        "full_compiler_flags":len(root_kernels),
        "primitive_inverse_clocks_per_flag":2,
        "A3_frames":len(a3_count),
        "native_compiler_C12_character":c_character,
        "cubic_torus_candidates":len(sigma_candidates),
        "matching_cubic_points":len(matching),
        "matching_carrier_Weyl_orbits":1,
        "sigma_torus_integral_residue_mod3":sigma_representative,
        "sigma_torus_standard_coordinates":[str(t) for t in sigma_torus],
        "analytic_consequence":"The whole cyclic native compiler is unitary-inner-conjugate to ct=kappa*exp(2pi i sigma_torus(0)); ct^9=kappa and ct^4=sigma_torus. The native-J fixed VOA is a conjugate of V_D5 tensor V_A3.",
        "proof_only":["compact E8 torus conjugacy and no diagram outer automorphisms",
                      "E8 lattice VOA is generated by its weight-one affine E8 algebra",
                      "existence of an inner conjugator, not a selected microscopic map"],
        "scope":"Full unitary E8 field automorphisms and transported Cartan/compiler charts. Conjugator need not preserve the original lattice Cartan, marking or sign-only lift subgroup.",
        "fundamental_physical_source":"NOT DERIVED",
        "source_sha256":old["PINS"],
        "promotion":False,
    }
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,default=HERE/"native-clock-conjugacy-results.json")
    args = parser.parse_args()
    result = calculate()
    args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
