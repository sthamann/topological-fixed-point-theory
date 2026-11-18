"""Synthesize the observable positive transfer from ordered return moments."""
import sympy as sp

from language import positive
from reconstruction import matrix, quotient, require


def block_moment(moments,level,shift):
    return sp.BlockMatrix([[moments[i+j+shift] for j in range(level)] for i in range(level)]).as_explicit()


def rows(m):
    return [[str(sp.simplify(v)) for v in row] for row in m.tolist()]


def reconstruct_return(moments,holdout=()):
    require(len(moments)>=4 and len(moments)%2==0,"supply an even number of training moments, starting at k=0")
    size = moments[0].rows
    require(size>0 and all(m.shape==(size,size) and m==m.H for m in [*moments,*holdout]),
            "all ordered moments must use one Hermitian marked basis")
    n = len(moments)//2
    gram = block_moment(moments,n,0)
    stepped = block_moment(moments,n,1)
    if not all(positive(q) for q in (gram,stepped,gram-stepped)):
        return {"status":"REFUTED_POSITIVE_CONTRACTION",
                "tested_forms":["Hankel", "shifted Hankel", "Hankel minus shifted Hankel"]}
    rank = gram.rank()
    if rank==0:
        if any(any(m) for m in holdout):
            return {"status":"HELD_OUT_MOMENT_DISTINGUISHES","first_difference":next(
                len(moments)+i for i,m in enumerate(holdout) if any(m))}
        return {"status":"ZERO_ACCESSIBLE_RETURN","rank":0}
    lower = block_moment(moments,n-1,0).rank()
    if lower!=rank:
        return {"status":"WORD_SPACE_STILL_GROWING","ranks":[lower,rank],
                "next_moments":[len(moments),len(moments)+1]}
    try:
        pivots,metric,ops = quotient(gram,{"time_step":stepped})
    except ValueError as error:
        return {"status":"REFUTED_MOMENT_CLOSURE","reason":str(error)}
    step = ops["time_step"]
    start = metric.inv()*gram.extract(pivots,range(size))
    coefficients = start
    for k,expected in enumerate([*moments,*holdout]):
        actual = (start.H*metric*coefficients).applyfunc(sp.expand)
        if actual!=expected:
            return {"status":"HELD_OUT_MOMENT_DISTINGUISHES" if k>=len(moments) else "REFUTED_MOMENT_CLOSURE",
                    "first_difference":k,"prediction":rows(actual),"observed":rows(expected)}
        coefficients = step*coefficients
    return {"status":"EXACT_ACCESSIBLE_POSITIVE_TRANSFER","rank":rank,
            "metric":rows(metric),"transfer":rows(step),"insertion":rows(start),
            "step_characteristic_polynomial":str(step.charpoly().as_expr()),
            "matched_training_moments":len(moments),"matched_held_out_moments":len(holdout),
            "predicted_next_moment":rows(start.H*metric*coefficients),
            "scope":"unique up to metric unitary on the cyclic span of the supplied return marks"}


def from_input(data):
    require(data.get("input_role")=="declared_finite_oracle","declare the return observations as finite controls")
    return {"input_role":data["input_role"],
            **reconstruct_return([matrix(v) for v in data["moments"]],
                                 [matrix(v) for v in data.get("held_out",[])])}
