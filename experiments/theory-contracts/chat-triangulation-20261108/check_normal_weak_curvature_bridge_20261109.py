"""Original normal clutch and weak-orbit connection, checked on the R032 cut."""
import contextlib
import io
import json
import runpy
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
import sympy as sp

repo = Path(sys.argv[1]).resolve()
directory = repo / "experiments/theory-contracts/chat-triangulation-20261108"
previous = directory / "check_original_whole_two_strips_20261109.py"
argv = sys.argv
try:
    sys.argv = [str(previous), str(repo)]
    with contextlib.redirect_stdout(io.StringIO()):
        existing = runpy.run_path(str(previous))
finally:
    sys.argv = argv

# Exact north/south degree-one monopole, with the actual clutch exp(i phi).
theta = sp.symbols("theta",real=True)
r = sp.symbols("r",positive=True)
north = sp.I*(1-sp.cos(theta))/2
south = -sp.I*(1+sp.cos(theta))/2
assert sp.simplify(north-south) == sp.I
curvature = sp.diff(north,theta)
assert sp.integrate(2*sp.pi*curvature,(theta,0,sp.pi)) == 2*sp.pi*sp.I
area = 4*sp.pi*r**2
projective_density = sp.simplify(curvature/(2*r**2*sp.sin(theta)))
assert projective_density == sp.I/(4*r**2)
curvature_floor = sp.simplify(area*sp.Abs(projective_density)**2)
assert curvature_floor == sp.pi/(4*r**2)
assert sp.exp(sp.I*sp.pi*(1-sp.cos(sp.pi))/2) == -1
x,y = sp.symbols("x y")
fplus,fminus = 2*sp.pi*sp.I*(x-y),2*sp.pi*sp.I*y
cross = sp.expand(((fplus-fminus)/2)**2*2/(8*sp.pi**2)).coeff(x,1).coeff(y,1)
assert cross == 1
tau = sp.symbols("tau",real=True)
product_energy = sp.pi**2*(sp.exp(tau)+4*sp.exp(-tau))
assert sp.simplify(sp.diff(product_energy,tau).subs(tau,sp.log(2))) == 0
assert sp.simplify(sp.diff(product_energy,tau,2).subs(tau,sp.log(2))) == 4*sp.pi**2

native = existing["native"]
source = existing["source"]
history = [np.asarray(h) for h in existing["history"]]
durations = existing["durations"]
identity = existing["identity"]
two = np.eye(2)
swap = source["swap"]
e = native.IDX[("e","R",0,0)]
nu = native.IDX[("nu","R",0,0)]


def norm(x):
    return float(np.linalg.norm(x,2))


def weak_parity(a):
    q = native.rep(1j,-1j*a,1j*np.eye(3))
    return q@swap@q.conj()@swap


projectors = [np.outer(h,h.conj())/np.vdot(h,h).real for h in history]
marks = [weak_parity(2*p-two) for p in projectors]
assert all(abs(q[e,e]+1)<1e-13 and abs(q[nu,nu]-1)<1e-13 for q in marks)
links = []
phase_links = []
link_comparison = []
for i,(p,q) in enumerate(zip(projectors,projectors[1:])):
    overlap = np.sqrt(np.trace(p@q).real)
    # Direct rotation = actual split Berry/Kato transport along the shortest
    # ray geodesic. This interpolation is an explicit comparison premise.
    u = (q@p+(two-q)@(two-p))/overlap
    assert norm(u.conj().T@u-two)<1e-13
    assert abs(np.linalg.det(u)-1)<1e-13
    assert norm(u@p@u.conj().T-q)<1e-13
    distance = np.sqrt(max(0,2-2*overlap))
    assert abs(norm(u-two)-distance)<1e-12
    link = source["real_action"](u)
    assert norm(marks[i+1]@link-link@marks[i])<1e-13
    links.append(link)
    left, right = [h / np.linalg.norm(h) for h in history[i:i+2]]
    q_left, q_right = [np.column_stack(([h[1].conjugate(),-h[0].conjugate()],h))
                       for h in (left,right)]
    phase_link = q_right @ q_left.conj().T
    assert np.linalg.norm(phase_link@left-right)<1e-13
    assert norm(phase_link@p@phase_link.conj().T-q)<1e-13
    phase_lift = source["real_action"](phase_link)
    assert norm(marks[i+1]@phase_lift-phase_lift@marks[i])<1e-13
    phase_links.append(phase_lift)
    overlaps = [np.vdot(right,v@left) for v in (two,u,phase_link)]
    assert abs(abs(overlaps[1])-1)<1e-13 and abs(overlaps[2]-1)<1e-13
    link_comparison.append({"cut_time":sum(durations[:i+1]),
                            "minimum_change_from_original_identity_weak_link":float(distance),
                            "ray_transport_phase":float(np.angle(np.vdot(right,left))),
                            "ray_transport_vector_defect":float(np.linalg.norm(u@left-right)),
                            "vector_transport_distance_from_original_link":norm(phase_link-two),
                            "native_ray_vector_complex_overlaps":[[float(z.real),float(z.imag)]
                                                                  for z in overlaps]})
for path in (links,phase_links):
    assert max(norm(a-b.conj().T) for a,b in zip(path,path[::-1]))<1e-13


def half(m, chosen_links, branch="A"):
    k=identity.copy()
    sea=0.
    for i,dt in enumerate(durations[:3]):
        h,_,_=existing["generator"](i,branch,m)
        assert norm(marks[i]@h-h@marks[i])<1e-13
        ev,u=eigh(h)
        if i:
            k=chosen_links[i-1]@k
        k=((u*np.exp(-dt*ev))@u.conj().T)@k
        sea+=dt*ev[:len(ev)//2].sum()
    return k,sea


def determinant(x):
    sign,logabs=np.linalg.slogdet(x)
    return sign*np.exp(logabs)


def occupations(b):
    x=identity+b
    idx = np.arange(len(x))
    z=determinant(x)
    def minor(removed):
        keep=idx[~np.isin(idx,removed)]
        return determinant(x[np.ix_(keep,keep)])
    he,hn,hen=minor([e]),minor([nu]),minor([e,nu])
    return np.array([z,z-he,z-hn,he+hn-2*hen])


def read(b):
    return (occupations(b)+np.array([1,-1,1,1])*occupations(b@marks[0]))/2


def response(m, chosen_links):
    k,sea=half(m,chosen_links)
    b=k.conj().T@k
    stamped=k.conj().T@marks[2]@k
    reference,ref_sea=existing["direct_whole"]("S",m)
    ref_half,ref_half_sea=half(m,chosen_links,"S")
    assert norm(ref_half.conj().T@ref_half-reference)<1e-12
    assert abs(2*ref_half_sea-ref_sea)<1e-12
    reference_value=read(reference)[0]
    values=np.array([read(b),read(stamped)])
    assert np.max(np.abs(values.imag)) < 1e-10*np.max(np.abs(values.real))
    return ((values*np.exp(2*sea-ref_sea)/reference_value).real,
            norm(marks[2]@k-k@marks[0]))


data={}
for m in (-.001,.001):
    native_values,native_defect=response(m,[identity]*len(links))
    orbit_values,orbit_defect=response(m,links)
    phase_values,phase_defect=response(m,phase_links)
    native_ratio=native_values[1]/native_values[0]
    orbit_ratio=orbit_values[1]/orbit_values[0]
    assert np.max(np.abs(orbit_ratio-[1,-1,1,1]))<1e-12
    assert np.max(np.abs(phase_values[1]/phase_values[0]-orbit_ratio))<1e-12
    assert abs(native_ratio[1]+1)>.5
    assert orbit_defect<1e-11
    data[str(m)]={"normal_sphere_floor":float(curvature_floor.subs(r,np.exp(m))),
                  "original_links_half_charge_intertwining_defect":native_defect,
                  "orbit_links_half_charge_intertwining_defect":orbit_defect,
                  "original_links_gauge_marked_over_plain":native_ratio.tolist(),
                  "orbit_links_gauge_marked_over_plain":orbit_ratio.tolist(),
                  "orbit_links_relative_empty_and_words_plain":orbit_values[0].tolist(),
                  "orbit_links_relative_empty_and_words_marked":orbit_values[1].tolist(),
                  "vector_links_half_charge_intertwining_defect":phase_defect,
                  "vector_links_gauge_marked_over_plain":(phase_values[1]/phase_values[0]).tolist(),
                  "vector_links_relative_empty_and_words_plain":phase_values[0].tolist(),
                  "vector_links_plain_over_ray_links":(phase_values[0]/orbit_values[0]).tolist()}
print(json.dumps({
    "scope":"geometric normal-sphere theorem plus conditional orbit transport on R032; not full original metric or QTheta",
    "normal_sphere_sharp_floor_for_area_A":"pi^2/A in -tr(X^2)/2 norm",
    "normalized_projective_holonomy_south_endpoint":"-I2",
    "mixed_normal_orbit_Chern_pairing":int(cross),
    "product_metric_minimum_area_ratio_W_over_N":2,
    "product_metric_minimum_energy_and_log_ratio_hessian":"4*pi^2",
    "links":link_comparison,
    "readout_order":["empty","right_e","right_nu","neutral_density_difference_squared"],
    "metrics":data,
},indent=2))
