"""Independent exact finite Gaussian source jets, including mixed contacts."""
import sympy as sp

from reconstruction import matrix, require, scalar
from source import MARKS


def gaussian_jet(data):
    require(data.get("input_role")=="declared_finite_oracle","Gaussian jet input requires a finite oracle")
    matrices = [matrix(data[name]) for name in ("kinetic","U","Y","UY")]
    n = matrices[0].rows
    require(1<=n<=6 and all(a.shape==(n,n) for a in matrices),"use four square matrices of equal size up to six")
    base,u,y,uy = matrices
    weight = list(map(scalar,data["scalar_weight_jet"]))
    reference = scalar(data["reference"])
    require(len(weight)==4 and reference!=0,"need four scalar weight jets and a fixed nonzero reference")

    s,t = sp.symbols("u y")
    determinant = (base+s*u+t*y+s*t*uy).det(method="berkowitz")
    ordered = [sp.expand(determinant).coeff(s,i).coeff(t,j) for i,j in ((0,0),(1,0),(0,1),(1,1))]
    direct = [base.det(),0,0,0]
    contact = sp.S.Zero
    for i in range(n):
        for mark,operator in enumerate((u,y,uy),1):
            block = base.copy()
            block[:,i] = operator[:,i]
            entry = block.det()
            direct[mark] += entry
            if mark==3:
                contact += entry
        for j in range(n):
            if i!=j:
                block = base.copy()
                block[:,i],block[:,j] = u[:,i],y[:,j]
                direct[3] += block.det()
    require(all(sp.expand(a-b)==0 for a,b in zip(ordered,direct)),"determinant jet and column expansion disagree")
    d0,du,dy,duy = direct
    w0,wu,wy,wuy = weight
    raw = list(map(sp.simplify,(w0*d0,w0*du+wu*d0,w0*dy+wy*d0,w0*duy+wu*dy+wy*du+wuy*d0)))
    z0,zu,zy,zuy = raw
    connected = sp.cancel(zuy/z0-zu*zy/z0**2) if z0 else None
    return {"status":"EXACT_FINITE_GAUSSIAN_JET","input_role":data["input_role"],
            "unnormalized_numerators":dict(zip(MARKS,map(str,raw))),
            "values":dict(zip(MARKS,(str(sp.simplify(v/reference)) for v in raw))),
            "kinetic_mixed_contact_numerator":str(sp.simplify(w0*contact)),
            "connected_UY":str(connected) if connected is not None else None,
            "base_kinetic_singular":d0==0,"physical_source_derived":False,
            "scope":"finite Gaussian determinant with supplied scalar-weight jet; no full R113 ensemble"}
