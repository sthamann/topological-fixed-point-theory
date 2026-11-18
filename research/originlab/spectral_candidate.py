"""A declared finite spectral-action candidate and its first source tests."""
from itertools import product
import json

from flint import acb_mat, arb, arb_mat, ctx

from r113 import ball
from reconstruction import require, scalar
from source import MARKS, NATIVE, source_refs
from verify import digest


def adjoint(m):
    return acb_mat([[m[j,i].conjugate() for j in range(m.nrows())] for i in range(m.ncols())])


def join(a,b,c,d):
    n = a.nrows()
    return acb_mat([[block[i,j] for block in pair for j in range(n)]
                    for pair in ((a,b),(c,d)) for i in range(n)])


def real(value):
    if not value.imag.contains(0) or not value.real.is_finite():
        raise ArithmeticError("expected a finite certified real value")
    return value.real


def half_blocks(data):
    from flint import acb
    y = {name:acb_mat([[acb(ball(v[0]),ball(v[1])) for v in row] for row in data[name]])
         for name in ("Yu","Yd","Ye","Ynu","MR")}
    zero = acb_mat(3,3)
    blocks = [(acb_mat(3,3),y[name],mult) for name,mult in (("Yu",6),("Yd",6),("Ye",2))]
    blocks.append((join(zero,zero,adjoint(y["MR"]),zero),
                   join(y["Ynu"],zero,zero,y["Ynu"].transpose()),1))
    for m,v,mult in blocks:
        n = m.nrows()
        z = acb_mat(n,n)
        diag = lambda sign:acb_mat([[sign if i==j else 0 for j in range(n)] for i in range(n)])
        yield join(z,m,adjoint(m),z),join(z,v,adjoint(v),z),join(diag(1),z,z,diag(-1)),mult


def strict_positive_gram(rows):
    # Interval Cholesky in LDL form. Exact singularity is reported as
    # unresolved except when a negative diagonal or minor is certified.
    if any(rows[i][i]<0 for i in range(len(rows))):
        i = next(i for i in range(len(rows)) if rows[i][i]<0)
        return {"status":"REFUTED_SCOPED","leading_order":i+1,"negative_value":str(rows[i][i])}
    for i in range(len(rows)):
        for j in range(i):
            if rows[i][i]*rows[j][j]-rows[i][j]**2<0:
                return {"status":"REFUTED_SCOPED","minor_indices":[j,i],
                        "negative_value":str(rows[i][i]*rows[j][j]-rows[i][j]**2)}
    if all(v==rows[0][0] for row in rows for v in row):
        return {"status":"POSITIVE_DEFINITE"} if len(rows)==1 and rows[0][0]>0 else {
            "status":"POSITIVE_SEMIDEFINITE","rank":int(rows[0][0]>0)}
    work = [row.copy() for row in rows]
    for i in range(len(work)):
        pivot = work[i][i]
        if pivot<0:
            return {"status":"REFUTED_SCOPED","leading_order":i+1,"negative_pivot":str(pivot)}
        if not pivot>0:
            return {"status":"UNRESOLVED","leading_order":i+1,"pivot":str(pivot)}
        for j in range(i+1,len(work)):
            for k in range(i+1,len(work)):
                work[j][k] -= work[j][i]*work[i][k]/pivot
    return {"status":"POSITIVE_DEFINITE"}


def negative_direction(rows,certificate):
    if "negative_pivot" not in certificate:
        return certificate
    size = certificate["leading_order"]
    block = arb_mat([row[:size-1] for row in rows[:size-1]])
    column = arb_mat([[-row[size-1]] for row in rows[:size-1]])
    try:
        solution = block.solve(column)
    except ZeroDivisionError:
        return certificate
    for bits in (12,24,48,96,144):
        scale = 2**bits
        rounded = [(solution[i,0].mid()*scale+arb("1/2")).floor().unique_fmpz()
                   for i in range(size-1)]
        if None in rounded:
            continue
        vector = [arb(a)/scale for a in rounded]+[arb(1)]
        form = sum(vector[i]*rows[i][j]*vector[j] for i,j in product(range(size),repeat=2))
        if form<0:
            certificate.update(witness_indices=list(range(size)),
                               witness_coefficients=[str(scalar(f"{a}/{scale}")) for a in rounded]+["1"],
                               witness_quadratic_form=str(form))
            break
    return certificate


def finite_candidate(root, data):
    """Sum all half-history pairs with a declared reference and ordered mark."""
    require(data.get("input_role")=="declared_alternative",
            "spectral candidate must be explicitly declared as an alternative")
    fields,masses = ([scalar(v) for v in data[key]] for key in ("half_higgs","half_masses"))
    require(1<=len(fields)<=6 and len(fields)==len(masses) and len(set(fields))==len(fields) and
            all(x.is_Rational for x in fields+masses) and all(m>0 for m in masses),
            "declare one to six distinct radial half fields and positive base masses")
    moment,width,cutoff,coupling = (scalar(data[key]) for key in
                                   ("momentum","half_time","cutoff","spectral_strength"))
    require(all(x.is_Rational for x in (moment,width,cutoff,coupling)) and
            moment>0 and width>0 and cutoff>0 and coupling>=0,
            "positive momentum, time and cutoff, nonnegative spectral strength required")
    path = root/NATIVE
    fixture = json.loads(path.read_text())
    bindings = source_refs(root,fixture["source_originals"])
    with ctx.workprec(192):
        p,w,lam,g = (ball(str(x)) for x in (moment,width,cutoff,coupling))
        xs = [ball(str(x)) for x in fields]
        base = [ball(str(x)) for x in masses]
        # Cache every half transfer. The 2x2 block exponential supplies the
        # exact Frechet derivative for the late Y insertion.
        cached = []
        for m,v,grading,multiplicity in half_blocks(fixture["couplings"]):
            n = grading.nrows()
            zero = acb_mat(n,n)
            identity = acb_mat([[int(i==j) for j in range(n)] for i in range(n)])
            y = -1j*grading*v
            rows = []
            for x in [*xs,arb(0)]:
                d = m+v*x
                h = p*grading-1j*grading*d
                exp_derivative = join(-w*h,-w*x*y,zero,-w*h).exp()
                rows.append((d,acb_mat([[exp_derivative[i,j] for j in range(n)] for i in range(n)]),
                    acb_mat([[exp_derivative[i,j+n] for j in range(n)] for i in range(n)])))
            cached.append((p*grading,identity,multiplicity,rows))
        require(sum(r[1].nrows()*r[2] for r in cached)==96,"native block multiplicities changed")

        def entry(i,j):
            trace,fermion,mark = arb(0),arb(0),arb(0)
            for bond,identity,multiplicity,rows in cached:
                dx,x,_ = rows[i]
                dy,y,dy_transfer = rows[j]
                operator = join(dx,bond,bond,dy)
                trace += multiplicity*real((-(operator*operator)/(lam**2)).exp().trace())
                core = identity+x*y
                determinant = real(core.det())
                if not determinant>0:
                    raise ArithmeticError("fermion determinant positivity could not be certified")
                fermion += multiplicity*determinant.log()
                mark -= multiplicity*real((core.inv()*x*dy_transfer).trace())/w
            return trace,fermion,mark

        ref_s,ref_f,_ = entry(-1,-1)
        size = len(fields)
        marks = [arb(0) for _ in MARKS]
        bare = [[arb(0) for _ in xs] for _ in xs]
        joint = [[arb(0) for _ in xs] for _ in xs]
        histories = []
        for i,j in product(range(size),repeat=2):
            s,f,y = entry(i,j)
            action = -g*(s-ref_s)
            weight = (action+f-ref_f).exp()
            if i<=j:
                bare[i][j]=bare[j][i]=action.exp()
                joint[i][j]=joint[j][i]=weight
            elif not weight.overlaps(joint[j][i]):
                raise ArithmeticError("reflection symmetry failed numerical enclosure")
            mass = base[i]*base[j]*weight
            u = xs[i]**2
            for k,term in enumerate((mass,mass*u,mass*y,mass*u*y)):
                marks[k] += term
            histories.append({"early":str(fields[i]),"late":str(fields[j]),"spectral_action":str(g*s),
                              "log_fermion":str(f),"late_y":str(y),"relative_kernel":str(weight)})
        connected = marks[3]/marks[0]-marks[1]*marks[2]/(marks[0]**2)
        return {"status":"CERTIFIED_DECLARED_FINITE_SOURCE",
                "input_role":"declared_alternative","dimension":96,"two_slice_dirac_dimension":192,
                "native_fixture_sha256":digest(path),"native_fixture_sources":bindings,
                "coordinate_space":"ordered pairs of supplied real scalar Higgs amplitudes",
                "operator":"D(x,y)=[[D_F(x),p*Gamma_F],[p*Gamma_F,D_F(y)]]",
                "measure":"positive product of declared finite half masses, multiplied by exp(-g*Tr exp(-D^2/Lambda^2))",
                "fermion":"ordinary full96 exterior determinant det(1+exp(-s h_x)exp(-s h_y))",
                "reference":"same empty x=y=0 weighted determinant",
                "parameter_inputs":data,"values":dict(zip(MARKS,map(str,marks))),
                "connected_UY":str(connected),"connected_excludes_zero":not connected.contains(0),
                "reflection":{"spectral_weight_only":negative_direction(bare,strict_positive_gram(bare)),
                              "combined_kernel":negative_direction(joint,strict_positive_gram(joint))},
                "history_receipts":histories,
                "unbounded_lebesgue_normalization":{"status":"REFUTED_SCOPED",
                    "spectral_bound":"bounded nonnegative f gives exp(-Tr f) >= exp(-dim(D)*sup f)",
                    "fermion_bound":"XY is similar to the positive matrix sqrt(X) Y sqrt(X), so det(1+XY) >= 1",
                    "scope":"bare bounded spectral weight or its product with this ordinary two-slice determinant on a flat infinite-volume base"},
                "physical_original_values":dict.fromkeys(MARKS),
                "scope":"two-slice radial surrogate, fixed diagnostic Yukawas and Majorana, ordinary determinant, no physical selector or geometry limit"}


def default_inputs():
    return {"input_role":"declared_alternative","half_higgs":[0,"7/20"],"half_masses":["1/2","1/2"],
            "momentum":"1/4","half_time":"2/5","cutoff":1,"spectral_strength":"1/4"}


def hostile_inputs():
    return {**default_inputs(),"half_higgs":[0,"7/20","7/10","21/20","7/5"],
            "half_masses":["1/5"]*5,"half_time":"1/200","spectral_strength":"1/100"}


def spectral_controls(root):
    small = finite_candidate(root,default_inputs())
    larger = hostile_inputs()
    results = {name:finite_candidate(root,{**larger,"spectral_strength":strength})
               for name,strength in (("unweighted",0),("reflection_counterexample","1/100"),
                                     ("stronger_profile","1/10"))}
    return {"status":"FINITE_SPECTRAL_CONTROLS","small":small,**results}
