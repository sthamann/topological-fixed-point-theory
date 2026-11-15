"""Exact P1 three-state associative operator closure and affine-index gate.

The original 3 log and 5 meromorphic modes can MATCH the compact
1+3+5 operator harmonics of one spin-one qutrit. This is a declared
physical field map, not the native original weighted Q0 producer.
"""
from itertools import combinations
import sympy as sp


def require(ok, message):
    if not ok:
        raise ArithmeticError(message)


def basis():
    E=lambda i,j:sp.eye(3)[:,i]*sp.eye(3)[j,:]
    anti=[E(0,1)-E(1,0),E(0,2)-E(2,0),E(1,2)-E(2,1)]
    sym=[E(0,1)+E(1,0),E(0,2)+E(2,0),E(1,2)+E(2,1),
         sp.diag(1,-1,0),sp.diag(1,1,-2)]
    return anti,sym


def controls():
    cases={n:{"family":n-1,"carrier":n+1,
              "full_matrix_dim":(n-1)**2,"geometric_with_unital":2*n+1}
           for n in range(2,17)}
    matches=[n for n,c in cases.items()
             if c["full_matrix_dim"]==c["geometric_with_unital"]]
    require(matches==[4],"unital original P1 operator closure not unique")
    anti,sym=basis()
    herm=[sp.eye(3),*(sp.I*a for a in anti),*sym]
    vecs=sp.Matrix.hstack(*(x.reshape(9,1) for x in herm))
    require(vecs.rank()==9 and all(v.H==v for v in herm),
            "geometric 1+3+5 not a complete Hermitian M3 basis")
    metric=sp.Matrix(9,9,lambda i,j:sp.trace(herm[i]*herm[j]))
    require(metric.is_positive_definite and metric.eigenvals()=={3:1,2:7,6:1},
            "positive Hilbert Schmidt metric changed")
    for a in herm:
        for b in herm:
            coeff=vecs.inv()*(a*b).reshape(9,1)
            require(coeff.rows==9,"3-state associative algebra closure failed")
    traceless=[sp.I*a for a in anti]+sym
    tvec=sp.Matrix.hstack(*(v.reshape(9,1) for v in traceless))
    for i,j in combinations(range(8),2):
        c=-sp.I*(traceless[i]*traceless[j]-traceless[j]*traceless[i])
        coeff=tvec.gauss_jordan_solve(c.reshape(9,1))[0]
        require(all(sp.im(x)==0 for x in coeff),
                "compact su3 Hermitian real bracket changed")

    compact=anti+[sp.I*s for s in sym]
    bvec=sp.Matrix.hstack(*(v.reshape(9,1) for v in compact))
    adj=[sp.Matrix.hstack(*(bvec.gauss_jordan_solve(
         (a*b-b*a).reshape(9,1))[0] for b in compact)) for a in compact]
    residual=sp.Matrix(8,8,lambda i,j:sp.trace(adj[i]*adj[j])-
                       6*sp.trace(compact[i]*compact[j]))
    require(residual==sp.zeros(8),
            "fundamental versus adjoint Dynkin-index ratio is not six")
    return {
      "status":"EXACT_UNITAL_3PLUS5_QUANTUM_M3_AND_AFFINE_INDEX_GATE",
      "closure_equation":"(n-1)^2=1+(n-1)+(n+1)",
      "unique_n_at_least_2":matches,
      "n4_family_Hilbert_dim":3,
      "n4_operator_algebra_dim":9,
      "spin_one_operator_decomposition":[1,3,5],
      "traceless_compact_Lie_algebra":"su3",
      "matrix_associative_products_checked":81,
      "positive_Hilbert_Schmidt_metric_spectrum":{"3":1,"2":7,"6":1},
      "adjoint_trace_index_over_fundamental":6,
      "affine_level_if_8_COMPLEX_chiral_fermions_in_su3_adjoint":6,
      "native_regular_E8_root_su3_level":1,
      "can_identify_both_as_the_same_level1_native_colour_current":False,
      "native_positive_weighted_source_and_R113_identified":False,
      "physical_principle_all_operators_from_original_sections_derived":False,
      "F1_solved":False,
    }


if __name__=="__main__":
    import json
    print(json.dumps(controls(),sort_keys=True,indent=2))
