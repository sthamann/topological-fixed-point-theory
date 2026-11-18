"""Compare actual C58 occupied-block det2, local subtractions and radial modes."""

import contextlib
import io
import json
from pathlib import Path
import runpy
import sys

import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
from scipy.special import zeta


repo = Path(sys.argv[1]).resolve()
base = repo / "experiments/theory-contracts/chat-triangulation-20261108"
with contextlib.redirect_stdout(io.StringIO()):
    radial = runpy.run_path(str(base / "check_tfpt_c58_relative_radial_measure_20261110.py"))
old = radial["old"]
ell = radial["ell"]
gamma = old["gamma"]
m, v = old["clifford_majorana"], old["clifford_higgs"]
baseline_schemes = json.loads((base / "tfpt-c58-renorm-results-20261110.json").read_text())["prescriptions"]


def log_cap_block(p, mass2):
    energy = np.hypot(p, np.sqrt(mass2))
    return (ell*mass2/(energy+p)
            + np.log1p(-mass2/(2*energy*(energy+p))*(-np.expm1(-2*ell*energy))))


def occupied_square(radius, sign):
    occupied = sign * np.diag(gamma) < 0
    leg = (m+radius*v)[np.ix_(occupied,~occupied)]
    return leg@leg.conj().T


references = {}
for sign, shift in ((1,.25),(-1,.75)):
    ref = occupied_square(0.,sign)
    values,u = np.linalg.eigh(ref)
    values = values.clip(0)
    # Majorana zero eigenspace plus the three actual nonzero eigenvalues.
    group = [values < 1e-12] + [np.arange(len(values))==j for j in np.flatnonzero(values>=1e-12)]
    references[sign] = (shift,ref,u,group,np.array([values[g].mean() for g in group]))


def logdet2(radius, cutoff):
    answer = 0.
    for sign,(shift,reference,u0,groups,reference_masses) in references.items():
        square = occupied_square(radius,sign)
        eigenvalues,u = np.linalg.eigh(square)
        overlaps = abs(u.conj().T@u0)**2
        weights = np.stack([overlaps[:,group].sum(axis=1) for group in groups],axis=1)
        p = (np.arange(cutoff)+shift)[:,None]
        a = log_cap_block(p,eigenvalues.clip(0))
        b = log_cap_block(p,reference_masses)
        d = a[:,:,None]-b[:,None,:]
        answer += np.sum((d-np.expm1(d))*weights)
        # Remaining O(p^-3) terms verified by cutoff doubling.
        answer -= ell**2/8 * np.linalg.norm(square-reference)**2 * zeta(2,cutoff+shift)
    return float(answer)


def direct_one_mode(radius, p):
    sign = 1 if p>0 else -1
    reference = occupied_square(0.,sign)
    values,u = np.linalg.eigh(reference)
    k = (u*np.exp(-log_cap_block(abs(p),values.clip(0))/2))@u.conj().T
    square = occupied_square(radius,sign)
    values,u = np.linalg.eigh(square)
    a = (u*np.exp(log_cap_block(abs(p),values.clip(0))))@u.conj().T
    matrix = k@a@k
    return np.linalg.slogdet(matrix)[1] - np.trace(matrix-np.eye(len(matrix))).real


def cartesian_tangent_upper():
    # log a_p(x) is concave by the positive Robin/Dirichlet determinant
    # representation. Subtract its linear ultraviolet term before summing.
    value = 0.
    for sign,(shift,_,u0,groups,reference_masses) in references.items():
        p = (np.arange(512)+shift)[:,None]
        x = reference_masses
        energy = np.hypot(p,np.sqrt(x))
        vacancy = x/(2*energy*(energy+p))
        heat = np.exp(-2*ell*energy)
        derivative = (-ell*x/(2*p*energy*(energy+p))
                      +(p/(4*energy**3)*(heat-1)-vacancy*ell/energy*heat)
                      /(1+vacancy*(heat-1)))
        occupied = sign*np.diag(gamma)<0
        higgs_leg = v[np.ix_(occupied,~occupied)]
        quadratic = np.sum(abs(u0.conj().T@higgs_leg)**2,axis=1)
        value += np.sum(derivative*np.array([quadratic[g].sum() for g in groups]))
    return float(value)


def det2_radial_tail_bound(radius):
    lower = []
    for sign,shift in ((1,.25),(-1,.75)):
        occupied = sign*np.diag(gamma)<0
        for index, charge in old["blocks"]:
            if charge==0:
                continue
            leg = v[np.ix_(index[occupied[index]],index[~occupied[index]])]
            for singular in np.linalg.svd(leg,compute_uv=False):
                for p in np.arange(12)+shift:
                    a = .5*np.exp(ell*(singular*radius-p))
                    if a>1:
                        lower.append((a-1-np.log(a),ell*singular*(a-1)))
    constant = sum(a for a,_ in lower)
    slope = sum(b for _,b in lower)
    # r^3 volume factor with a linear lower tangent of a-1-log(a).
    tail = np.log(2*np.pi**2 * (radius**3/slope + 3*radius**2/slope**2
                              + 6*radius/slope**3 + 6/slope**4))-constant
    return float(tail)


for r in (.1,.55):
    expected = direct_one_mode(r,.25) + direct_one_mode(r,-.75)
    tail = sum(ell**2/8*np.linalg.norm(occupied_square(r,sign)-ref)**2*zeta(2,1+shift)
               for sign,(shift,ref,*_) in references.items())
    assert abs(logdet2(r,1)+tail-expected)<1e-12


def density_and_peak(log_weight):
    def density(r):
        return 2*np.pi**2*r**3*np.exp(log_weight(r))
    mass, _ = quad(density,0,8,epsabs=1e-11,points=[.15,.4,1.])
    sample = np.linspace(.001,2,81)
    i = int(np.argmax([density(r) for r in sample]))
    assert 0<i<len(sample)-1
    peak = minimize_scalar(lambda r:-3*np.log(r)-log_weight(r),
                           bounds=(sample[i-1],sample[i+1]),method="bounded",
                           options={"xatol":1e-9})
    return {"mass_0_8":mass,"normalized_cartesian_null_density":1/mass,
            "radial_shell_peak_r":peak.x,
            "log_cartesian_weight_at_shell_peak":log_weight(peak.x)}


rows = {}
for name, record in baseline_schemes.items():
    if name.startswith("C59"):
        continue
    c = record["c"]
    log_weight = lambda r:radial["relative_log"](r,512)-c*r*r
    rows[name] = {"c":c,**density_and_peak(log_weight)}
    assert abs(rows[name]["normalized_cartesian_null_density"]
               -record["normalized_cartesian_null_density"])<1e-5

x = np.array([v["c"] for v in rows.values()])
y = np.log([v["normalized_cartesian_null_density"] for v in rows.values()])
order = np.argsort(x)
slopes = np.diff(y[order])/np.diff(x[order])
assert np.all(np.diff(slopes)<0)
log_det2 = lambda r:logdet2(r,512)
chosen = density_and_peak(log_det2)
assert abs(chosen["normalized_cartesian_null_density"] -
           density_and_peak(lambda r:logdet2(r,1024))["normalized_cartesian_null_density"]) < .005
samples = {}
for r in (.05,.1,.2,.35,.55):
    v2 = logdet2(r,2048)
    coarse = logdet2(r,1024)
    old_weight = radial["relative_log"](r,1024)
    samples[str(r)] = {"modewise":old_weight,"det2":v2,"det2_cutoff_difference":abs(coarse-v2),
                       "finite_difference_divided_by_r_squared":(old_weight-v2)/r**2}
j = {str(step):2*logdet2(step,1024)/step**2 for step in (.004,.002,.001)}
zero_hessian_c = baseline_schemes["modewise_subtracted"]["null_log_weight_radial_hessian"]/2
cartesian_upper = cartesian_tangent_upper()
assert cartesian_upper < min(x)
log_tail_bound = det2_radial_tail_bound(8.)
assert log_tail_bound < -30., log_tail_bound
print(json.dumps({"old_schemes":rows,
                  "null_log_density_sorted_secant_slopes":slopes.tolist(),
                  "maximum_absolute_affine_fit_residual":float(np.max(abs(y-np.polyval(np.polyfit(x,y,1),x)))),
                  "modewise_c_required_to_remove_r_squared_term":zero_hessian_c,
                  "modewise_log_cartesian_weight_tangent_upper_per_r_squared":cartesian_upper,
                  "det2_log_tail_mass_upper_at_r_8":log_tail_bound,
                  "det2_Majorana_reference":{**chosen,"local_hessian_checks":j,"r_samples":samples},
                  "scope":"Old C58 fixed H/Majorana quarter cap, det2 of normalized occupied operator relative to H=0 Majorana, flat d4H. No primitive R113 amplitude."},indent=2))
