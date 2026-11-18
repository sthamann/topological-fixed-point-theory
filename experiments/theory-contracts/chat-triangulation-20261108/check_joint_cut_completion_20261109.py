"""Complete original R032 diagnostic halfword Gram blocks by selector sector."""
import contextlib
import io
import json
import runpy
import sys
from pathlib import Path

import numpy as np

repo = Path(sys.argv[1]).resolve()
base = repo / "experiments/theory-contracts/chat-triangulation-20261108"
prior = base / "check_normal_weak_curvature_bridge_20261109.py"
argv = sys.argv
try:
    sys.argv = [str(prior),str(repo)]
    with contextlib.redirect_stdout(io.StringIO()):
        old = runpy.run_path(str(prior))
finally:
    sys.argv = argv

src = old["source"]
v = old["native"]
identity, mid = old["identity"], old["marks"][2]
assert np.linalg.norm(mid.conj().T-mid)<1e-13
assert np.linalg.norm(mid@mid-identity)<1e-13
plus, minus = (identity+mid)/2, (identity-mid)/2
# Use the pre-existing formal weak sigma_x, with its part transverse
# to the actual midpoint Higgs ray. This is an insertion, not a new link.
px = np.array([[0.,1.],[1.,0.]])
p = old["projectors"][2]
w = p@px@(np.eye(2)-p)+(np.eye(2)-p)@px@p
a = v.rep(0,1j*w,np.zeros((3,3)))
swap = v.Kswap()
q = -1j*(a+swap@a.conj()@swap)
assert np.linalg.norm(mid@q@mid+q)<1e-13
weak = q @ src["grading"]
h0 = old["history"][0]
coupling = src["inner_form"](h0)
left = v.rep(0,np.eye(2),np.zeros((3,3)))
# Do not use left as a physical projection; the v252 Yukawa right columns
# already lie on the corresponding left particle legs.
assert np.linalg.norm(left@coupling[:,[old["e"],old["nu"]]]-coupling[:,[old["e"],old["nu"]]])<1e-13
start_e, start_nu = coupling[:,old["e"]], coupling[:,old["nu"]]


def norm(x):
    return float(np.linalg.norm(x,2))


def psd_factor(a):
    lam, u = np.linalg.eigh((a+a.conj().T)/2)
    assert min(lam)>-1e-13
    keep = lam>1e-13*max(1.,max(lam))
    return np.sqrt(lam[keep])[:,None]*u[:,keep].conj().T


def complex_pair(z):
    return [float(z.real),float(z.imag)]


rows = {}
for m in (-.001,.001):
    k,sea = old["half"](m,[identity]*len(old["links"]))
    k_ref,ref_sea = old["existing"]["direct_whole"]("S",m)
    norm_scalar = np.exp(2*sea-ref_sea)/old["read"](k_ref)[0].real
    fields = np.column_stack((k@start_e,weak@k@start_nu))
    full = fields.conj().T@fields
    raw = [fields.conj().T@s@fields for s in (plus,minus)]
    for mat in raw:
        a0,b0 = mat[0,0].real,mat[1,1].real
        re_probe,im_probe = np.array([1,1]),np.array([1,1j])
        re_read = np.vdot(re_probe,mat@re_probe).real
        im_read = np.vdot(im_probe,mat@im_probe).real
        polar = (re_read-a0-b0 + 1j*(a0+b0-im_read))/2
        assert abs(polar-mat[0,1])<1e-13
    assert norm(full-sum(raw))<1e-12
    a,b = np.linalg.eigh(full)
    assert a[0]>0
    h = (b/np.sqrt(a))@b.conj().T
    u = fields@h
    effect = h@raw[0]@h
    leakage = (identity-u@u.conj().T)@plus@u
    assert norm(effect-effect@effect-leakage.conj().T@leakage)<1e-12
    factors = [psd_factor(x) for x in raw]
    reconstruction = np.concatenate(factors,axis=0)
    assert norm(reconstruction.conj().T@reconstruction-full)<1e-12
    rank_total = len(a)
    ranks = [len(x) for x in factors]
    proj0 = np.diag([1]*ranks[0]+[0]*ranks[1])
    assert norm(reconstruction.conj().T@proj0@reconstruction-raw[0])<1e-12
    caps = [np.sqrt((a_[0,0]*a_[1,1]).real) for a_ in raw]
    c = full[0,1]
    angle = c/abs(c)
    lower = max(-caps[0],abs(c)-caps[1])
    upper = min(caps[0],abs(c)+caps[1])
    assert lower<upper
    for x in (lower,upper):
        z = x*angle
        assert abs(z)<=caps[0]+1e-12
        assert abs(c-z)<=caps[1]+1e-12
        # Hold the entire unheated Gram and both sector diagonals fixed.
        # These are admissible Gram completions, not new physical sources.
        p_alt = np.diag(np.diag(raw[0])).astype(complex)
        p_alt[0,1],p_alt[1,0]=z,z.conjugate()
        q_alt = full-p_alt
        assert np.max(np.abs(np.diag(q_alt)-np.diag(raw[1])))<1e-12
        assert min(np.linalg.eigvalsh(p_alt))>-1e-12
        assert min(np.linalg.eigvalsh(q_alt))>-1e-12
    rows[str(m)] = {
        "scope":"conditional original R032 one-particle halfword, candidate weak mark, unchanged native links; no genuine QTheta or Padm",
        "unheated_cross":complex_pair(c),"candidate_even_cross":complex_pair(raw[0][0,1]),
        "same_unheated_and_sector_diagonals_allow_cross_line_interval":[float(lower),float(upper)],
        "line_phase":complex_pair(angle),"sector_cross_radii":list(map(float,caps)),
        "original_halfword_span_rank":rank_total,
        "minimal_selector_closed_rank":sum(ranks),
        "sector_ranks":ranks,
        "compressed_even_effect_eigenvalues":list(map(float,np.linalg.eigvalsh(effect))),
        "defect_from_idempotence_operator_norm":norm(effect-effect@effect),
        "polarization_reconstructs_complex_cross_with_two_positive_diagonal_probes":True,
        "retained_even_electron_halfword_fraction":float((raw[0][0,0]/full[0,0]).real),
        "retained_even_transverse_neutrino_halfword_fraction":float((raw[0][1,1]/full[1,1]).real),
        "same_diagnostic_S_normalizer_for_all_candidate_marks":float(norm_scalar),
    }
print(json.dumps({"scope":"exact Gram completion and two actual original diagnostic metrics on unchanged R032; not full weighted source",
                  "sector_words":["Y_e on initial right electron","Y_nu on initial right neutrino followed by transverse field at cut"],
                  "rows":rows},indent=2))
