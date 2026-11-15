"""Recover and certify nonnegative masses for a declared family of native-form atoms."""
import sympy as sp
import z3

from reconstruction import require, scalar
from source import MARKS, check_fibre_pushforward, complex_pair, fibre_pushforward_data, iter_weighted_fibres, rational


def number(value):
    require(value.is_Rational, "inverse fibre coefficients must be exact rationals")
    return z3.RealVal(str(value))


def linear(row, variables):
    return z3.Sum([number(a)*x for a,x in zip(row,variables) if a])


def solution(solver, variables):
    verdict = solver.check()
    if verdict == z3.unsat:
        return None
    if verdict != z3.sat:
        raise ArithmeticError(f"rational solver returned {verdict}")
    model = solver.model()
    values = [model.eval(x,model_completion=True) for x in variables]
    return sp.Matrix(len(values),1,[sp.Rational(v.numerator_as_long(),v.denominator_as_long()) for v in values])


def equality_system(a,b):
    x = [z3.Real(f"mass_{i}") for i in range(a.cols)]
    solver = z3.SolverFor("QF_LRA")
    solver.add(*[v>=0 for v in x])
    solver.add(*[linear(a.row(i),x)==number(b[i]) for i in range(a.rows)])
    return solver,x


def certificate(a,b,support=None):
    """Find and independently replay a Farkas or maximal-support functional."""
    y = [z3.Real(f"dual_{i}") for i in range(a.rows)]
    singleton = support is None
    solver = z3.SolverFor("QF_LRA")
    solver.add(linear(b,y)==(-1 if singleton else 0))
    for i in range(a.cols):
        value = linear(a.col(i),y)
        solver.add(value>=0 if singleton else value==0 if i in support else value>=1)
    result = solution(solver,y)
    if result is None or (b.T*result)[0] != (-1 if singleton else 0):
        raise ArithmeticError("missing dual certificate")
    coefficients = a.T*result
    if not all(v>=0 if singleton else v==0 if i in support else v>=1
               for i,v in enumerate(coefficients)):
        raise ArithmeticError("invalid dual certificate")
    return {"dual":list(map(str,result)), "atom_values":list(map(str,coefficients)),
            "target_value":str((b.T*result)[0])}


def geometry(a,b):
    solver,x = equality_system(a,b)
    first = solution(solver,x)
    if first is None:
        return {"feasible":False,"certificate":certificate(a,b)},None,None
    total = first
    samples = 1
    while True:
        missing = [i for i,w in enumerate(total) if not w]
        if not missing:
            break
        solver.push()
        solver.add(z3.Or([x[i]>0 for i in missing]))
        witness = solution(solver,x)
        solver.pop()
        if witness is None:
            break
        total += witness
        samples += 1
    centre = total/samples
    support = [i for i,w in enumerate(centre) if w]
    if a*centre != b or any(w<0 for w in centre):
        raise ArithmeticError("invalid primal certificate")
    subspace = a[:,support]
    rank = subspace.rank()
    report = {"feasible":True,"masses":list(map(str,centre)), "support":support, "active_rank":rank,
              "mass_freedom":len(support)-rank}
    if len(support)<a.cols:
        report["support_certificate"] = certificate(a,b,set(support))
    return report,centre,subspace.nullspace()


def lower_bound(a,b,c):
    base,x = equality_system(a,b)
    solver = z3.Optimize()
    solver.add(*base.assertions())
    handle = solver.minimize(linear(c,x))
    point = solution(solver,x)
    if point is None:
        raise ArithmeticError("known feasible optimization became infeasible")
    if str(handle.lower()) == "-1*oo":
        ray_solver,ray_vars = equality_system(a.col_join(c),sp.zeros(a.rows,1).col_join(sp.Matrix([-1])))
        ray = solution(ray_solver,ray_vars)
        if ray is None or a*ray != sp.zeros(a.rows,1) or (c*ray)[0] != -1 or any(v<0 for v in ray):
            raise ArithmeticError("invalid unbounded ray")
        return {"value":"-inf","ray":list(map(str,ray))}
    value = (c*point)[0]
    y = [z3.Real(f"bound_{i}") for i in range(a.rows)]
    dual = z3.SolverFor("QF_LRA")
    dual.add(linear(b,y)==number(value))
    dual.add(*[linear(a.col(i),y)<=number(c[i]) for i in range(a.cols)])
    multiplier = solution(dual,y)
    if (multiplier is None or a*point != b or any(v<0 for v in point) or
            (b.T*multiplier)[0] != value or any(v<0 for v in c.T-a.T*multiplier)):
        raise ArithmeticError("invalid optimum certificate")
    return {"value":str(value),"masses":list(map(str,point)),"dual":list(map(str,multiplier))}


def interval(a,b,c):
    lower = lower_bound(a,b,c)
    upper = lower_bound(a,b,-c)
    upper["value"] = "inf" if upper["value"]=="-inf" else str(-sp.Rational(upper["value"]))
    return {"minimum":lower,"maximum":upper,"maximum_certificate_uses_negative_objective":True}


def assemble(level, ids):
    size = level["fine_size"]
    time = level["time"]
    u_time,y_time = map(rational,(time["t_U"],time["t_Y"]))
    reference = complex_pair(level["fine_reference_scalar"])
    require(reference != 0,"a bound nonzero absolute reference is required")
    scale = complex_pair(level["fine_scalar_factor"])/reference
    data = level["coarse_graining"]
    atoms = level["fine_atoms"]
    require(type(size) is int and size>0 and bool(time["id"]) and 0<=u_time<=y_time,
            "valid common word dimension and ordered time required")
    require(all("original_measure_mass" not in atom for atom in atoms), "fine masses must remain unknown")
    unit = [{**atom,"original_measure_mass":[1,0]} for atom in atoms]
    fine = dict(iter_weighted_fibres(unit,size,u_time,y_time))
    require(set(fine)==set(ids), "snapshots must transport exactly the same atom identities")
    embedding,coarse,assignment,coarse_scale = fibre_pushforward_data(
        data,ids,size,time["id"],u_time,y_time,scale)
    fit = set(level.get("fit_marks",["Z0","ZU","ZY"]))
    require(fit <= set(MARKS) and "Z0" in fit and fit != set(MARKS),
            "each snapshot must fit Z0 and retain at least one withheld mark")
    rows = {"fit":[],"withheld":[]}
    for cell,forms in coarse.items():
        for mark in MARKS:
            columns = [scale*embedding.H*fine[k][mark]*embedding if assignment[k]==cell
                       else sp.zeros(embedding.cols) for k in ids]
            target = coarse_scale*forms[mark]
            for i in range(target.rows):
                for j in range(target.cols):
                    parts = [column[i,j].as_real_imag() for column in columns]
                    for half,name in enumerate(("real","imaginary")):
                        coeff = [v[half] for v in parts]
                        value = target[i,j].as_real_imag()[half]
                        if any(coeff) or value:
                            rows["fit" if mark in fit else "withheld"].append(
                                ([time["id"],cell,mark,i,j,name],coeff,value))
    def replay(masses, support, kernel, equations, target):
        records = [{**atom,"original_measure_mass":[str(masses[ids.index(atom["original_configuration_id"])]),0]}
                   for atom in atoms]
        forms = dict(iter_weighted_fibres(records,size,u_time,y_time))
        checked = check_fibre_pushforward(data,forms,size,time["id"],u_time,y_time,scale)
        if checked["status"] != "EXACT_SUPPLIED_MARKED_PUSHFORWARD":
            raise ArithmeticError("reconstructed masses fail independent fibre replay")
        zero_form = sum((a["Z0"] for a in forms.values()),sp.zeros(size))
        coordinate = None
        for direction in kernel:
            delta = sum((w*fine[ids[i]]["Z0"] for i,w in zip(support,direction)),sp.zeros(size))
            entry = next(((i,j) for i in range(size) for j in range(size) if delta[i,j]),None)
            if entry is not None:
                i,j = entry
                half = 0 if sp.re(delta[i,j]) else 1
                coeff = sp.Matrix([[fine[name]["Z0"][i,j].as_real_imag()[half] for name in ids]])
                coordinate = {"entry":[i,j],"component":("real","imaginary")[half],
                              **interval(equations,target,coeff)}
                break
        return {"time_id":time["id"],"candidate_fibre_records":records,"replay":checked,
                "candidate_fine_Z0_before_outer_contraction":
                    [[[str(x) for x in value.as_real_imag()] for value in row] for row in zero_form.tolist()],
                "fine_Z0_fixed_by_full_constraints":coordinate is None,
                "fine_Z0_separating_coordinate":coordinate}
    return rows,replay


def run(data):
    levels = data["snapshots"]
    require(isinstance(levels,list) and levels and len({l["time"]["id"] for l in levels})==len(levels),
            "supply distinct time snapshots")
    require(len({s["coarse_graining"]["reference_binding"] for s in levels})==1,
            "time snapshots must share one absolute reference binding")
    def bindings(level):
        return {a["original_configuration_id"]:(a["whole_fiber_binding"],a["original_measure_binding"])
                for a in level["fine_atoms"]}
    bound = bindings(levels[0])
    require(bound and all(bindings(level)==bound for level in levels),
            "transport needs matching fibre and measure bindings at each time")
    ids = list(bound)
    fit,withheld,replays = [],[],[]
    for level in levels:
        rows,replay = assemble(level,ids)
        fit.extend(rows["fit"])
        withheld.extend(rows["withheld"])
        replays.append(replay)
    if "total_mass" in data:
        fit.append((["independent_absolute_total_mass"],[sp.S.One]*len(ids),rational(data["total_mass"])))

    def system(rows):
        return (sp.Matrix([row for _,row,_ in rows]) if rows else sp.zeros(0,len(ids)),
                sp.Matrix([value for _,_,value in rows]) if rows else sp.zeros(0,1))

    a,b = system(fit)
    report,centre,kernel = geometry(a,b)
    result = {"atom_order":ids,"fit_rows":[label for label,_,_ in fit],"fit":report,
              "physical_source_identity_proved":False}
    if centre is None:
        return {**result,"status":"FIT_EXCLUDED_IN_DECLARED_ATOM_FAMILY"}
    support = report["support"]
    probes = []
    for row in data.get("independent_probes",[]):
        require(set(row["coefficients"])<=set(ids), "unknown atom in independent probe")
        probes.append(([row["id"]],[scalar(row["coefficients"].get(i,0)) for i in ids],scalar(row["expected"])))
    varying = next((row for row in withheld if any(
        (sp.Matrix([row[1]])[:,support]*direction)[0] for direction in kernel)),None)
    result["withheld_forms_predicted_from_fit"] = varying is None
    if varying:
        result["separating_form"] = {"coordinate":varying[0],"observed":str(varying[2]),
                                     **interval(a,b,sp.Matrix([varying[1]]))}
    result["independent_probe_intervals"] = [
        {"coordinate":label,"observed":str(target),**interval(a,b,sp.Matrix([coeff]))}
        for label,coeff,target in probes]

    # Infer using fit first; replay the joint originals in a second, explicit feasibility test.
    validation_rows = fit+withheld+probes
    validation_a,validation_b = system(validation_rows)
    full,weights,full_kernel = geometry(validation_a,validation_b)
    result["validation_rows"] = [label for label,_,_ in validation_rows]
    result["validation"] = full
    if weights is None:
        return {**result,"status":"WITHHELD_EXCLUDED_IN_DECLARED_ATOM_FAMILY"}
    signatures = {}
    for i,name in enumerate(ids):
        column = validation_a.col(i)
        scale = next((abs(v) for v in column if v),sp.S.One)
        signatures.setdefault(tuple(column/scale),[]).append({"atom":name,"column_scale":str(scale)})
    result["proportional_atom_classes"] = [group for group in signatures.values() if len(group)>1]
    result["candidate_snapshots"] = [replay(weights,full["support"],full_kernel,validation_a,validation_b)
                                     for replay in replays]
    return {**result,"status":"UNIQUE_MASSES_IN_DECLARED_ATOM_FAMILY" if full["mass_freedom"]==0
            else "MASSES_UNDERDETERMINED_IN_DECLARED_ATOM_FAMILY"}
