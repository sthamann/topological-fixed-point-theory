"""Scoped R113 identifiability checks using the original R032 conditional kernels."""
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import json

from flint import acb, acb_mat, arb, ctx

from source import NATIVE

WINDOWS = Path("experiments/theory-contracts/chat-triangulation-20261108/"
               "TFPT-R113-GETRENNTE-ZEITFENSTER-20261109.txt")


def ball(value):
    # The native domain loader converts these JSON numbers to binary floats.
    return arb(str(Fraction(value)))


def squared_singular_values(rows):
    y = acb_mat([[acb(ball(z[0]), ball(z[1])) for z in row] for row in rows])
    adjoint = acb_mat([[y[j, i].conjugate() for j in range(y.nrows())]
                       for i in range(y.ncols())])
    values = (adjoint * y).eig(multiple=True)
    if len(values) != y.ncols() or any(not ev.imag.contains(0) or not ev.real > 0 for ev in values):
        raise ArithmeticError("could not certify a positive charged singular spectrum")
    return [ev.real for ev in values]


def pair_trace(squared, early, late, momentum, width):
    ea = (momentum**2 + early**2*squared).sqrt()
    eb = (momentum**2 + late**2*squared).sqrt()
    a = (width*ea).cosh()
    b = (width*ea).sinh()/ea
    c = (width*eb).cosh()
    d = (width*eb).sinh()/eb
    cross = momentum**2 + early*late*squared
    eb_derivative = late**2*squared/eb
    d_derivative = (width*c-d)*eb_derivative/eb
    trace = 2*(a*c+b*cross*d)
    derivative = 2*(a*width*eb*d*eb_derivative + b*(early*late*squared*d + cross*d_derivative))
    return trace, derivative, ea+eb


def conditional_word(spectra, early, late, momentum, width):
    # Each charged singular pair has a determinant-one transfer. At z=1 the
    # exterior factor is tr R+2. For the conjugate colour charges at z^3=1,
    # their product is (tr R-1)^2; neutral leptons still use z=1.
    weights, scores = [arb(1), arb(1)], [arb(0), arb(0)]
    sea = arb(0)
    for colour, squares in spectra:
        multiplicity = 6 if colour else 2
        for squared in squares:
            trace, derivative, energy = pair_trace(squared, early, late, momentum, width)
            sea -= width*multiplicity*energy
            for i, shift in enumerate((2, -1 if colour else 2)):
                factor = trace+shift
                if not factor > 0:
                    raise ArithmeticError("could not certify the charged exterior factor")
                weights[i] *= factor**multiplicity
                scores[i] += multiplicity*derivative/factor
    k0 = (weights[0]+2*weights[1])*sea.exp()/3
    mean_y = -(weights[0]*scores[0]+2*weights[1]*scores[1])/(width*(weights[0]+2*weights[1]))
    return k0, mean_y


def witness(root, early_radii=("7/20","2/5"), late="-11/20", momentum="1/4", time="4/5"):
    """Same whole radial law and three scalar marks, different mixed response.

    This fixes the conditional kernels on four given R032 histories. The
    two base measures are witnesses of ambiguity, not original path masses.
    """
    from reconstruction import scalar
    values = [scalar(v) for v in (*early_radii,late,momentum,time)]
    if len(values)!=5 or not all(x.is_Rational for x in values):
        raise ValueError("two rational radii, a late field, momentum and duration required")
    a,b,late,momentum,time = map(str,values)
    a,b = map(Fraction,(a,b))
    if not 0 < a < b or Fraction(late) == 0 or Fraction(momentum) <= 0 or Fraction(time) <= 0:
        raise ValueError("two increasing positive early radii and nonzero late field required")
    domain = root/NATIVE
    raw = domain.read_bytes()
    fixture = json.loads(raw)
    if fixture["domain_id"] != "original-variable-H-quarter-full96-diagnostic-v1":
        raise ValueError("unexpected diagnostic coupling dictionary")
    for ref in fixture["source_originals"]:
        if sha256((root/ref["path"]).read_bytes()).hexdigest()!=ref["sha256"]:
            raise ValueError(f"changed original input: {ref['path']}")
    with ctx.workprec(192):
        spectra = [(name!="Ye",squared_singular_values(fixture["couplings"][name]))
                   for name in ("Yu","Yd","Ye")]
        width = ball(time)/2
        p = ball(momentum)
        histories = (a,-a,b,-b)
        evaluated = [conditional_word(spectra,ball(x),ball(late),p,width) for x in histories]
        k = [row[0] for row in evaluated]
        y = [row[1] for row in evaluated]
        empty_ref,_ = conditional_word(spectra,arb(0),arb(0),p,width)
        delta_a,delta_b = y[0]-y[1],y[2]-y[3]
        if any(not x>0 for x in k) or delta_a.contains(0) or delta_b.contains(0):
            return {"status":"UNRESOLVED","reason":"orientation contrast or positivity not separated from zero"}
        u = [ball(x*x) for x in histories]
        # Exact effective-mass perturbation: each radial history has mass 1/2.
        direction = (delta_b,-delta_b,-delta_a,delta_a)
        epsilon = arb("1/1024")
        effective = [[arb("1/4")+sign*epsilon*d for d in direction] for sign in (1,-1)]
        if any(not value>0 for weights in effective for value in weights):
            return {"status":"UNRESOLVED","reason":"chosen witness step did not certify positive masses"}
        gap = 2*epsilon*(u[0]-u[2])*delta_a*delta_b
        if gap.contains(0):
            raise ArithmeticError("no separated mixed response")
        connected = [sum(m*x*v for m,x,v in zip(w,u,y))-
                     sum(m*x for m,x in zip(w,u))*sum(m*v for m,v in zip(w,y)) for w in effective]
        if not (connected[0]-connected[1]).overlaps(gap):
            raise ArithmeticError("independent direct contraction disagrees")
        normalized = [[sum(w),sum(m*x for m,x in zip(w,u)),sum(m*v for m,v in zip(w,y)),
                       sum(m*x*v for m,x,v in zip(w,u,y))] for w in effective]
        if any(not x.overlaps(v) for x,v in zip(normalized[0][:3],normalized[1][:3])):
            raise ArithmeticError("lower responses do not replay")
        return {"status":"EXCLUDED_SCOPED","source_id":"R032:charged84:fixed-connection:Z3",
                "native_domain_sha256":sha256(raw).hexdigest(),
                "original_reference":str(WINDOWS),
                "holds_fixed":["all four conditional exterior kernels","connection","local sea factor per history",
                               "one common outer scalar reference","whole weighted radial path law"],
                "matched_observations":["Z0","ZU","ZY","weighted_radial_history_law"],
                "separating_observation":"ZUY (equivalently connected UY)",
                "why":"Within each radial class the two masses sum to 1/2; "
                      "the signed Yukawa shifts cancel across classes.",
                "early_fields":list(map(str,histories)),"late_field":str(late),"time":str(time),
                "momentum":str(momentum),"late_y":[str(x) for x in y],
                "conditional_unmarked_kernels":[str(x) for x in k],
                "conditional_log_Z0_over_empty_reference":[str((x/empty_ref).log()) for x in k],
                "conditional_marks_over_common_reference":[
                    dict(zip(("Z0","ZU","ZY","ZUY"),[str(x) for x in row])) for row in normalized],
                "witness_step":"1/1024",
                "conditional_base_masses":[[str(v/x) for v,x in zip(w,k)] for w in effective],
                "conditional_connected":[str(x) for x in connected],
                "connected_plus_minus":[str(gap.lower()),str(gap.upper())],
                "formula":"2*epsilon*(u_A-u_B)*(y_Aplus-y_Aminus)*(y_Bplus-y_Bminus)",
                "missing_original_choice":"conditional charged orientation/history weights over the fixed radial path law",
                "scope":"conditional original R032 charged exterior diagnostic with its Z3 filter; "
                        "the two positive base measures are not the full R113 measure or selector",
                "physical_original_values":dict.fromkeys(("Z0","ZU","ZY","ZUY"))}
