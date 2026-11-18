"""Exact rectangular intertwiner constraints on supplied finite carriers."""
import sympy as sp

from language import positive
from reconstruction import matrix, require
from return_moments import rows


def solve(data):
    metrics = (matrix(data[key]) for key in ("domain_metric","target_metric"))
    relations = {label:tuple(matrix(pair[key]) for key in ("domain","target"))
                 for label,pair in data["relations"].items()}
    return solve_operators(*metrics,relations)


def solve_operators(left,right,relations):
    n,m = left.rows,right.rows
    require(0<n<=m and n*m<=4096 and all(positive(g) and g.det()!=0 for g in (left,right)),
            "need positive definite metrics and a bounded rectangular carrier")
    require(isinstance(relations,dict) and relations,"labelled operator relations required")
    adj_left,adj_right = left.inv(),right.inv()
    systems = []
    for label,(a,b) in relations.items():
        require(a.shape==(n,n) and b.shape==(m,m),"operator does not match its declared carrier")
        for suffix,x,y in (("",a,b),(" adjoint",adj_left*a.H*left,adj_right*b.H*right)):
            equation = sp.kronecker_product(sp.eye(m),x.T)-sp.kronecker_product(y,sp.eye(n))
            systems.append((label+suffix,equation))
    equations = sp.Matrix.vstack(*(eq for _,eq in systems))
    kernel = [vec.reshape(m,n) for vec in equations.nullspace()]
    result = {"linear_solution_dimension":len(kernel),"domain_dimension":n,"target_dimension":m,
              "physical_source_derived":False,"scope":"declared same-source finite adjoint-closed operator relations"}
    if not kernel:
        core = systems[:]
        for row in systems:
            reduced = [item for item in core if item is not row]
            if reduced and sp.Matrix.vstack(*(eq for _,eq in reduced)).rank()==n*m:
                core = reduced
        return {**result,"status":"EXCLUDED_RECTANGULAR_INTERTWINER","conflicting_relations":[label for label,_ in core]}
    attempts = kernel[:]
    attempts.extend(kernel[0]+k for k in kernel[1:])
    for v in attempts:
        if v.rank()==n:
            ratio = left.inv()*v.H*right*v
            # Positive metric adjunction implies V*V is in the domain
            # commutant. Polar normalization can be checked independently.
            spectral = ratio.eigenvals()
            if not all(value.is_positive for value in spectral):
                continue
            identity = sp.eye(n)
            inv_root = sp.zeros(n)
            for value in spectral:
                projector = identity
                for other in spectral:
                    if other!=value:
                        projector = projector*(ratio-other*identity)/(value-other)
                inv_root += projector/sp.sqrt(value)
            w = (v*inv_root).applyfunc(sp.simplify)
            require(verify_embedding(left,right,relations,w)["status"]=="EXACT_SUPPLIED_EMBEDDING",
                    "isometry replay failed")
            return {**result,"status":"EXACT_RECTANGULAR_ISOMETRY","map":rows(w)}
    return {**result,"status":"NO_ISOMETRY_IN_ONE_DIMENSIONAL_SPACE" if len(kernel)==1 and kernel[0].rank()<n
            else "UNRESOLVED_RECTANGULAR_ISOMETRY","basis":[rows(v) for v in kernel]}


def verify_embedding(left, right, operators, embedding):
    """Check a proposed same-source embedding in the actual quotient metrics."""
    require(embedding.shape==(right.rows,left.rows),"embedding has the wrong carrier dimensions")
    gram = (embedding.H*right*embedding-left).applyfunc(sp.simplify)
    defects = {"metric":gram.rank()}
    for name,(a,b) in operators.items():
        defects[name] = (embedding*a-b*embedding).applyfunc(sp.simplify).rank()
        defects[name+" adjoint"] = (embedding*left.inv()*a.H*left-right.inv()*b.H*right*embedding
                                   ).applyfunc(sp.simplify).rank()
    return {"status":"EXACT_SUPPLIED_EMBEDDING" if not any(defects.values()) else "REFUTED_SUPPLIED_EMBEDDING",
            "defect_ranks":defects,"physical_source_derived":False}
