"""Two-face source response, not an assumed local quantum instrument.

All rational values are read from the already specified quarter-circle
product collar, alpha=1/4 and ell=4 log 2. Infinite-shift formulas are checked
coefficientwise, without truncating the shift to a finite/cyclic matrix.
"""
import json
import sys
from pathlib import Path
import sympy as s


def ratio(n):
    return s.Rational(2)**(-(4*n+1))


def occupation_chart(n):
    return 1/(1+ratio(n)**2)


def graph(n):
    r=ratio(n)
    return s.Matrix([[1,r],[r,r*r]])/(1+r*r)


def exact_check():
    p=s.diag(1,0)
    rows=[]
    for n in range(-3,4):
        r=ratio(n); a=occupation_chart(n); g=graph(n)
        j=s.Matrix([1,r])/s.sqrt(1+r*r)
        null=s.Matrix([-r,1])/s.sqrt(1+r*r)
        assert g*g == g and g.T == g
        assert s.simplify(j*j.T-g) == s.zeros(2)
        assert g*null == s.zeros(2,1)
        leak=s.simplify((null.T*p*g*p*null)[0])
        assert leak == a*(1-a) and leak > 0
        m=graph(n+1)
        # diag(U,U) has identity 2x2 block from mode n to n+1.
        assert m != g  # it fails to commute with the fixed graph.
        block=s.simplify((s.Matrix([1,ratio(n+1)]).T /
                          s.sqrt(1+ratio(n+1)**2)*p*j)[0])
        assert s.simplify(block-s.sqrt(a*occupation_chart(n+1))) == 0
        forward=a*occupation_chart(n+1)
        backward=a*occupation_chart(n-1)
        comm=s.factor(forward-backward)
        assert comm > 0
        rows.append({'n':n,'a':str(a),'Lstar_L':str(forward),
                     'L_Lstar':str(backward),'commutator':str(comm),
                     'null_vector_face_mask_norm_squared':str(leak)})
    # Fixed-graph facewise unitary classification, linear intertwiner space:
    # the diagonal graph blocks force U0 to preserve each B eigenspace;
    # nonzero offdiagonal graph block then forces U1=U0. Check both simple
    # spectrum (complex dimension 3) and a double degeneracy (dimension 5).
    from sympy import linear_eq_to_matrix
    dimensions=[]
    for modes in [(-1,0,1),(0,0,1)]:
        aa=s.diag(*(occupation_chart(n) for n in modes))
        cc=s.diag(*(ratio(n)*occupation_chart(n) for n in modes))
        gg=aa.row_join(cc).col_join(cc.row_join(s.eye(3)-aa))
        variables=s.symbols('x:18')
        u0=s.Matrix(3,3,variables[:9]);u1=s.Matrix(3,3,variables[9:])
        ss=s.diag(u0,u1)
        coeff,_=linear_eq_to_matrix(list(ss*gg-gg*ss),variables)
        dim=18-coeff.rank();dimensions.append(dim)
        assert dim==(3 if len(set(modes))==3 else 5)
    a0=occupation_chart(0);a1=occupation_chart(1)
    w0=s.sqrt(a0*a1)
    # xi=(e0+e1)/sqrt2. This is a coherent boundary-vector comparison,
    # not a physically source-prepared quantum state claim.
    overlap=w0/2
    return {'verdict':'NONCOMMUTING_GEOMETRIC_RESPONSE_DERIVED_PHYSICAL_OPERATION_OPEN',
            'alpha':'1/4','ell':'4 log 2','mode_checks':rows,
            'coherent_overlap_I_L_on_two_mode_boundary_vector':str(overlap),
            'coherent_overlap_I_iL':str(s.I*overlap),
            'fixed_graph_intertwiner_dimensions_simple_degenerate':dimensions,
            'raw_face_mask_descends_to_G_quotient':False,
            'compressed_response_descends':True,
            'physical_local_instrument_derived':False,
            'physical_state_preparation_derived':False,
            'raw_P1_source_selected':False,
            'internal_carrier_or_coupling_assumed':False}


if __name__=='__main__':
    result=exact_check()
    target=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name('boundary_operation_response_results.json')
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
