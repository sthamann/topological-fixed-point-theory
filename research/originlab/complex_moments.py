"""Exact complex marked moment search on a declared discrete spectrum."""
import sympy as sp

from reconstruction import matrix, require, scalar
from return_moments import reconstruct_return


def analyze(data):
    levels = tuple(map(scalar,data["spectrum"]))
    n = len(levels)
    require(1<=n<=10 and len(set(levels))==n and all(x.is_Rational and 0<x<=1 for x in levels),
            "use one to ten distinct rational transfer values in (0,1]")
    diagonal = matrix(data["diagonal_observations"])
    require(diagonal.shape==(n,2) and all(x.is_Rational for x in diagonal),
            "provide n real moment rows with two marks")
    observations = data["mixed_observations"]
    require(isinstance(observations,dict) and all(k.isascii() and k.isdecimal() and len(k)<=2
                                                  for k in observations),"give mixed steps as decimal keys")
    mixed = {int(k):scalar(v) for k,v in observations.items()}
    require(len(mixed)==len(observations) and all(k<=64 and v.is_number and
            sp.re(v).is_Rational and sp.im(v).is_Rational for k,v in mixed.items()),
            "invalid mixed step or Gaussian rational observation")
    vandermonde = sp.Matrix([[value**k for value in levels] for k in range(n)])
    weights = vandermonde.inv()*diagonal
    if any(x<0 for x in weights):
        return {"status":"EXCLUDED_COMPLEX_MOMENTS","reason":"negative diagonal spectral weight"}
    caps_squared = [weights[j,0]*weights[j,1] for j in range(n)]
    active = [j for j,value in enumerate(caps_squared) if value>0]
    steps = list(mixed)
    step_rows = sp.Matrix(len(mixed),n,lambda i,j: levels[j]**steps[i])
    equations = step_rows.extract(range(len(mixed)),active)
    target = sp.Matrix(len(mixed),1,list(mixed.values()))
    seed, parameters = equations.gauss_jordan_solve(target) if (
        equations.rank()==equations.row_join(target).rank()) else (None,None)
    if seed is None:
        return {"status":"EXCLUDED_COMPLEX_MOMENTS","reason":"inconsistent mixed linear constraints"}
    if "feasible_seed" in data:
        supplied = sp.Matrix(list(map(scalar,data["feasible_seed"])))
        require(supplied.shape==(n,1) and step_rows*supplied==target and
                all(sp.re(z).is_Rational and sp.im(z).is_Rational for z in supplied) and
                all(sp.expand(z*sp.conjugate(z))<=cap for z,cap in zip(supplied,caps_squared)),
                "feasible seed must satisfy all moment and positivity constraints")
        seed = supplied.extract(active,[0])
    else:
        seed = seed.subs(dict.fromkeys(parameters,0))
    free = equations.nullspace()
    seed_feasible = all(sp.expand(z*sp.conjugate(z))<=caps_squared[j] for j,z in zip(active,seed))
    if not seed_feasible:
        return {"status":"UNRESOLVED_COMPLEX_FEASIBILITY" if free else "EXCLUDED_COMPLEX_MOMENTS",
                "reason":"particular solution lies outside a spectral positivity disk",
                "linear_freedom":len(free)}

    def expand(values):
        result = sp.zeros(n,1)
        for j,value in zip(active,values):
            result[j] = value
        return result

    result = {"status":"MIXED_MOMENTS_FIXED","linear_freedom":len(free),
              "mixed_spectral_weights":list(map(str,expand(seed))),
              "scope":"fixed discrete spectrum, at most two visible dimensions per level, no physical source claim"}
    if not free:
        return result
    # Seek a rational direction that fixes the saturated disks. If this
    # conservative search fails, boundary curvature may still permit motion.
    boundary = [i for i,j in enumerate(active) if sp.expand(seed[i]*sp.conjugate(seed[i]))==caps_squared[j]]
    fixed = sp.Matrix(len(boundary),len(active),lambda i,j:int(j==boundary[i]))
    directions = [expand(v) for v in equations.col_join(fixed).nullspace()]
    if not directions:
        return {**result,"status":"UNRESOLVED_COMPLEX_BOUNDARY_FREEDOM"}
    answer = next((k, direction) for k in range(n) for direction in directions
                  if (vandermonde[k,:]*direction)[0]!=0)
    k,direction = answer
    # A small rational step in the imaginary direction exposes phases
    # invisible to real-only readouts and preserves all known mixed data.
    delta = sp.I*direction
    base = expand(seed)
    step = next((sp.Rational(1,2**i) for i in range(33) if all(
        all(sp.expand(z*sp.conjugate(z))<=cap for z,cap in zip(base+s*delta/2**i,caps_squared))
        for s in (-1,1))),None)
    if step is None:
        return {**result,"status":"UNRESOLVED_COMPLEX_STEP_BOUND"}
    witnesses = []
    training = 2*(n+1)
    for sign in (-1,1):
        coeff = base+sign*step*delta
        moments = []
        for power in range(training+2):
            row = sp.Matrix([[value**power for value in levels]])
            a,b = list(row*weights)
            mark = sp.expand((row*coeff)[0])
            moments.append(sp.Matrix([[a,sp.conjugate(mark)],[mark,b]]))
        replay = reconstruct_return(moments[:training],moments[training:])
        require(replay["status"]=="EXACT_ACCESSIBLE_POSITIVE_TRANSFER","complex witness failed Hankel replay")
        witnesses.append({"mixed_spectral_weights":list(map(str,coeff)),"mixed_value":str(moments[k][1,0]),
                          "transfer_replay":replay})
    return {**result,"status":"MIXED_MOMENTS_UNDERDETERMINED",
            "next_probe":{"step":k,"component":"imaginary","word":["Y_adjoint",f"T^{k}","U"]},
            "witnesses":witnesses}


def control():
    return analyze({"spectrum":["1/2","1/4","1/8"],
                    "diagonal_observations":[[3,3],["7/8","7/8"],["21/64","21/64"]],
                    "mixed_observations":{"0":0}})
