"""Recompute frozen TFPT benchmarks and keep comparison schemes explicit."""
from decimal import Decimal, localcontext
from pathlib import Path
import json

from flint import arb, ctx

from verify import digest

REGISTRY = Path("verification/predictions_frozen.json")
ALPHA_PROOF = Path("verification/v382_alpha_quillen_exact.py")
ALPHA_MEASUREMENTS = (
    ("CODATA 2022","137.035999177","0.000000021","https://physics.nist.gov/cuu/pdf/JPCRD2022CODATA.pdf"),
    ("rubidium recoil","137.035999206","0.000000011","https://doi.org/10.1038/s41586-020-2964-7"),
    ("electron magnetic moment plus SM","137.035999166","0.000000015","https://doi.org/10.1103/PhysRevLett.130.071801"),
    ("caesium recoil","137.035999046","0.000000027","https://doi.org/10.1126/science.aap7706"),
)


def alpha_value(order, factor=2, exponent="-5/4"):
    """Interval bisection for the stated root and each nested ablation."""
    cc = 1/(8*arb.pi())
    base = 1/(6*arb.pi())
    cutoff = 48*cc**4

    def f(a):
        q = cutoff if order>=2 else arb(0)
        if order>=3:
            q *= (-factor*a).exp()
        if order>=4:
            q *= (arb(exponent)*(1-q).log()).exp()
        return a**3-(2*cc**3*a*a if order>=1 else arb(0))+(arb(164)/5)*cc**6*(base+q).log()

    lo,hi = arb("1/256"),arb("1/64")
    if not f(lo)<0 or not f(hi)>0:
        raise ArithmeticError("electromagnetic root lacks the expected bracket")
    for _ in range(122):
        mid = (lo+hi)/2
        value = f(mid)
        if value<0:
            lo = mid
        elif value>0:
            hi = mid
        else:
            break
    return 1/lo.union(hi)


def audit(root):
    registry = root/REGISTRY
    data = json.loads(registry.read_text())
    rows = {v["id"]:v for v in data["predictions"]}
    with ctx.workprec(192):
        codata = arb(ALPHA_MEASUREMENTS[0][1])
        sigma = arb(ALPHA_MEASUREMENTS[0][2])
        names = ("leading","quadratic","puncture","e_minus_2a","full")
        results = []
        for order,name in enumerate(names):
            value = alpha_value(order)
            results.append({"stage":name,"alpha_inverse":str(value),
                            "relative_to_codata":str((value-codata)/codata)})
        full = alpha_value(4)
        reference = arb(rows["ALPHA_INV"]["frozen_value"])
        if not abs(full-reference)<arb("1e-21"):
            raise ArithmeticError("computed value no longer matches frozen alpha")
        ablations = []
        for f,e in ((1,"-5/4"),(3,"-5/4"),(2,"-1"),(2,"-3/2")):
            val = alpha_value(4,f,e)
            ablations.append({"factor":f,"exponent":e,"alpha_inverse":str(val),
                              "signed_codata_sigma_measurement_only":str((val-codata)/sigma)})
        experimental = [{"measurement":name,"value":v,"sigma":s,"source":link,
                         "signed_sigma_measurement_only":str((full-arb(v))/arb(s))}
                        for name,v,s,link in ALPHA_MEASUREMENTS]

    # PDG 2026 pole masses in MeV. Propagate independent listed errors
    # for this diagnostic; the frozen ratios have a different scheme.
    with localcontext() as decimal_context:
        decimal_context.prec = 40
        electron,muon,tau = map(Decimal,("0.51099895069","105.6583755","1776.93"))
        se,sm,st = map(Decimal,("0.00000000016","0.0000023","0.09"))
        e_z,mu_z,tau_z = map(Decimal,("2.77713e-6","5.85042e-4","0.99378e-2"))
        leptons = []
        for key,measured,relative_sigma,running in (
            ("MMU_OVER_MTAU",muon/tau,((sm/muon)**2+(st/tau)**2).sqrt(),mu_z/tau_z),
            ("ME_OVER_MMU",electron/muon,((se/electron)**2+(sm/muon)**2).sqrt(),e_z/mu_z),
        ):
            row = rows[key]
            value = Decimal(row["frozen_value"])
            leptons.append({"id":key,"layer":row["layer"],"source_ratio":str(value),
                           "pole_ratio":str(measured),"relative_excess":str(value/measured-1),
                           "running_sm_at_MZ":str(running),"relative_excess_at_MZ":str(value/running-1),
                           "relative_measurement_sigma":str(relative_sigma),
                           "required_pole_over_source_transfer":str(measured/value),
                           "pole_identity":"REFUTED_SCOPED",
                           "scheme_matched_prediction":"BLOCKED_MISSING_SOURCE_TO_POLE_MAP"})
        cross_ratio = Decimal(rows["MMU_OVER_MTAU"]["frozen_value"])**2/Decimal(rows["ME_OVER_MMU"]["frozen_value"])
        measured_cross = (muon/tau)**2/(electron/muon)
    return {"status":"EMPIRICAL_COMPARISON_SCOPED","registry_sha256":digest(registry),
            "alpha":{"stages":results,"altered_profile":ablations,"measurements":experimental,
                     "equation_derived_from_original_source":False,
                     "open_target":"ALPHA.QUILLEN.EXACT.01","ledger_witness_sha256":digest(root/ALPHA_PROOF)},
            "leptons":{"pole_data":"https://pdg.lbl.gov/2026/listings/particle_properties.html",
                       "running_data":"https://doi.org/10.1103/fdcc-ycph",
                       "rows":leptons,"watchdog_relative_window":"0.03",
                       "source_seed_independent_ratio":str(cross_ratio),
                       "pole_seed_independent_ratio":str(measured_cross),
                       "same_seed_pole_map":"REFUTED_SCOPED"},
            "physical_gate_G7":"OPEN"}
