"""Exact identifiability of connected R113 marks from positive histories."""

from itertools import combinations
from math import comb

import sympy as sp

from reconstruction import matrix, require, scalar
from source import MARKS


def connected_history(data):
    histories = data["histories"]
    require(isinstance(histories,list) and histories and all(isinstance(label,str) and label for label in histories)
            and len(histories)==len(set(histories)), "distinct named histories required")
    n = len(histories)
    coefficients = data["mark_coefficients"]
    require(set(coefficients)==set(MARKS), "supply the four common unnormalized mark coefficients")
    marks = matrix([coefficients[key] for key in MARKS])
    require(marks.shape==(4,n) and all(v.is_Rational for v in marks),
            "use real rational coefficients with one column per history")
    empty = marks[0,:]
    require(all(v>0 for v in empty), "each declared history needs a strictly positive selected empty contraction")

    # q = p/(empty*p) makes relative observations linear and preserves
    # feasible distributions under the inverse map p=q/sum(q).
    constraints = data.get("relative_constraints",{})
    rows = [list(empty)]
    for item in constraints.values():
        row = matrix([item["coefficients"]])
        value = scalar(item["value"])
        require(row.cols==n and value.is_Rational and all(v.is_Rational for v in row),
                "use real rational relative observations")
        rows.append(list(row-value*empty))
    system = sp.Matrix(rows).row_join(sp.Matrix([1,*[0]*len(constraints)]))
    reduced, pivots = system.rref()
    if pivots[-1]==n:
        return {"status":"EXCLUDED_IN_DECLARED_RELATIVE_HISTORY_CLASS","reason":"inconsistent observations"}
    rank = len(pivots)
    if comb(n,rank)>4096:
        return {"status":"BOUND_EXCEEDED","basis_count":comb(n,rank)}
    rhs = reduced[:rank,n]
    vertices = set()
    for support in combinations(range(n),rank):
        basis = reduced[:rank,list(support)]
        if basis.det()==0:
            continue
        positive_part = basis.inv()*rhs
        if any(value<0 for value in positive_part):
            continue
        weights = [sp.S.Zero]*n
        for i,value in zip(support,positive_part):
            weights[i] = value
        vertices.add(tuple(weights))
    if not vertices:
        return {"status":"EXCLUDED_IN_DECLARED_RELATIVE_HISTORY_CLASS","reason":"no positive-rule completion"}

    def observables(weights):
        z,u,y,uy = marks*sp.Matrix(weights)
        require(z==1,"a normalized history completion must have empty response one")
        return u,y,uy,uy-u*y

    def witness(weights, values):
        total = sum(weights)
        return {"probabilities":dict(zip(histories,(str(v/total) for v in weights))),
                "relative_marks":dict(zip(MARKS,("1",*map(str,values[:3])))),
                "connected_UY":str(values[3])}

    data_rows = [(weights,observables(weights)) for weights in sorted(vertices)]
    u_values,y_values = ({row[1][i] for row in data_rows} for i in (0,1))
    constant = data_rows[0][1][3]
    report = {"absolute_reference_used":False,"vertices":len(vertices),
              "fixed_single_marks":{"U":len(u_values)==1,"Y":len(y_values)==1},
              "retained_relative_observations":list(constraints),
              "scope":"convex normalized completions of the declared finite positive history menu"}

    def next_query():
        if len(u_values)==1:
            u = next(iter(u_values))
            return {"word":["U-minus-mean","Y"],"fixed_mean":str(u),
                    "coefficients":list(map(str,marks[3,:]-u*marks[2,:])),
                    "normalization":"same_empty"}
        if len(y_values)==1:
            y = next(iter(y_values))
            return {"word":["U","Y-minus-mean"],"fixed_mean":str(y),
                    "coefficients":list(map(str,marks[3,:]-y*marks[1,:])),
                    "normalization":"same_empty"}
        return {"word":["U"],"coefficients":list(map(str,marks[1,:])),
                "normalization":"same_empty","then":"center U at the observed normalized mean and read its ordered product with Y"}

    first = data_rows[0]
    different = next((row for row in data_rows if row[1][3]!=constant),None)
    if different is not None:
        return {**report,"status":"CONNECTED_RESPONSE_UNDERDETERMINED",
                "next_observation":next_query(),"witnesses":[witness(*first),witness(*different)]}
    if len(u_values)==1 or len(y_values)==1:
        return {**report,"status":"CONNECTED_RESPONSE_FIXED","connected_UY":str(constant)}

    # The connected response is affine on the whole convex set only if
    # at least one of the single marks is constant. Check the midpoint
    # explicitly, since vertex values can otherwise all coincide.
    pair = next((row for row in data_rows if first[1][0]!=row[1][0] and first[1][1]!=row[1][1]),None)
    if pair is None:
        base = first[1]
        first,pair = (next(row for row in data_rows if row[1][i]!=base[i]) for i in (0,1))
    midpoint = tuple((a+b)/2 for a,b in zip(first[0],pair[0]))
    middle = observables(midpoint)
    offset = (first[1][0]-pair[1][0])*(first[1][1]-pair[1][1])/4
    require(middle[3]==constant+offset and offset!=0,"midpoint must break the equal vertex response")
    return {**report,"status":"CONNECTED_RESPONSE_UNDERDETERMINED",
            "all_vertex_connected_values":str(constant),"midpoint_offset":str(offset),
            "next_observation":next_query(),
            "witnesses":[witness(*first),witness(*pair),witness(midpoint,middle)]}


def control():
    return connected_history({"histories":["a","b"],
                             "mark_coefficients":{"Z0":[1,2],"ZU":[0,4],"ZY":[0,6],"ZUY":[0,12]}})
