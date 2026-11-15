"""Find exact separating moments of the native gauge-invariant Higgs norm."""
from collections import deque
from itertools import product
import json

import sympy as sp

from reconstruction import require
from source import NATIVE, rational, source_refs
from verify import digest


def norm_coefficients(proposal, genes):
    require(len(genes)==len(proposal["slots"]), "one gene per original Fourier slot required")
    allowed = list(map(rational,proposal["allowed_multipliers"]))
    require(all(type(i) is int and 0<=i<len(allowed) for i in genes), "invalid Fourier gene")
    field = {0:sp.Matrix(list(map(rational,proposal["base_constant"])))}
    zero = sp.zeros(len(field[0]),1)
    for gene,slot in zip(genes,proposal["slots"]):
        harmonic = slot["harmonic"]
        component = slot["component"]
        amplitude = allowed[gene]*rational(slot["unit"])
        if harmonic == 0:
            field[0][component] += amplitude
        else:
            factor = 1 if slot["kind"]=="cos" else -sp.I
            field.setdefault(harmonic,zero.copy())[component] += amplitude*factor/2
            field.setdefault(-harmonic,zero.copy())[component] += amplitude*sp.conjugate(factor)/2
    norm = {}
    for a,x in field.items():
        for b,y in field.items():
            norm[a+b] = norm.get(a+b,sp.S.Zero)+(x.T*y)[0]
    return {k:sp.expand(v) for k,v in norm.items() if k>=0}


def basis_rows(signatures, labels):
    """Close the constant under multiplication in the finite evaluation algebra."""
    size = len(signatures)
    generators = [tuple(row[j] for row in signatures) for j in range(len(labels))]
    queue = deque([((),(sp.S.One,)*size)])
    seen = {()}
    echelon = {}
    selected = []
    while queue and len(selected)<size:
        word,values = queue.popleft()
        remainder = list(values)
        for pivot,row in sorted(echelon.items()):
            if coefficient := remainder[pivot]:
                remainder = [a-coefficient*b for a,b in zip(remainder,row)]
        pivot = next((i for i,v in enumerate(remainder) if v),None)
        if pivot is None:
            continue
        echelon[pivot] = [v/remainder[pivot] for v in remainder]
        selected.append((word,values))
        for j,gen in enumerate(generators):
            child = tuple(sorted((*word,j)))
            if child not in seen:
                seen.add(child)
                queue.append((child,tuple(a*b for a,b in zip(values,gen))))
    require(len(selected)==size, "distinct Higgs-norm signatures must span their finite function algebra")
    return selected


def probe_plan(proposal):
    modes = range(len(proposal["allowed_multipliers"]))
    genes = list(product(modes,repeat=len(proposal["slots"])))
    values = [norm_coefficients(proposal,g) for g in genes]
    harmonics = sorted(set().union(*(set(v) for v in values)))
    labels = [(k,part) for k in harmonics for part in (0,1)]
    signatures = [tuple(v.get(k,sp.S.Zero).as_real_imag()[part] for k,part in labels) for v in values]
    indices = [j for j in range(len(labels)) if len({v[j] for v in signatures})>1]
    labels = [labels[j] for j in indices]
    groups = {}
    for gene,signature in zip(genes,signatures):
        key = tuple(signature[j] for j in indices)
        groups.setdefault(key,[]).append(list(gene))
    basis = basis_rows(list(groups),labels)
    evaluation = sp.Matrix([values for _,values in basis])
    return {
        "history_count":len(genes),
        "norm_signature_classes":list(groups.values()),
        "independent_moment_count_including_mass":len(basis),
        "measurement_plan":[{
            "word":[{"harmonic":labels[j][0],"component":("real","imaginary")[labels[j][1]]} for j in word],
            "value_per_class":list(map(str,values)),
        } for word,values in basis],
        "class_totals_from_ordered_measured_moments":
            [list(map(str,row)) for row in evaluation.inv().tolist()],
        "mass_reconstructible_from_these_moments":"one total per norm-signature class",
        "physical_source_derived":False,
    }


def native_probe_plan(root, domain):
    require(domain.resolve()==(root/NATIVE).resolve(),"the probe plan needs the original Q domain")
    data = json.loads(domain.read_text())
    refs = source_refs(root,data["source_originals"])
    proposal = data["field_proposal"]
    require(len(proposal["base_constant"])==4,"the original quaternionic Higgs chart has four real components")
    # Replay the original invisible harmonic refinement. It is a test
    # direction on a finer cutoff, not an assumption about measure support.
    reflected = {**proposal,"slots":[*proposal["slots"],
                {"kind":"cos","harmonic":8,"component":0,"unit":proposal["slots"][0]["unit"]}]}
    selector = []
    for j in range(len(proposal["allowed_multipliers"])):
        coeff = norm_coefficients(reflected,[2,1,1,j])
        selector.append({"gene":j,"U_7":str(coeff[7]),"U_16":str(coeff[16])})
    require(sp.Matrix([[1,sp.Rational(row["U_7"]),sp.Rational(row["U_16"])] for row in selector]).det()!=0,
            "original three-point gauge-invariant selector lost rank")
    return {"domain_sha256":digest(domain),"source_refs":refs,"native_proposal":probe_plan(proposal),
            "fine_refinement":[-7,6],"fine_harmonic":8,"fine_proposal":probe_plan(reflected),
            "three_point_selector":selector,
            "scope":"exact spatial field probes; weighted moments need the original common fibre and reference"}
