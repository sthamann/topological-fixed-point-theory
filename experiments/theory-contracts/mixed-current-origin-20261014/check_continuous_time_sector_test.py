#!/usr/bin/env python3
"""Finite controls for the written compactness theorem; no new sector claim."""
import hashlib
import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
import sympy as s

OUT = Path(__file__).resolve().parent
REPO = Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')


def need(ok, text):
    if not bool(ok):
        raise RuntimeError(text)


z = s.symbols('z', real=True)
levels = [0, 1, 2, 3, 8]
weights = [s.Rational(1, 8), s.Rational(1, 4), s.Rational(1, 4),
           s.Rational(1, 8), s.Rational(1, 4)]
laplace = sum(w * z ** e for w, e in zip(weights, levels))
tail = sum(w for w, e in zip(weights, levels) if e > 2)
margin = s.cancel((1 - laplace - (1-z**2)*tail) / (1-z))
need(all(c >= 0 for c in s.Poly(margin, z).all_coeffs()),
     'spectral tail margin must be positive throughout0<z<1')
need(sum(weights) == 1, 'normalized control state')
real_window=.4
real_average_defect=sum(float(w)*(1-math.sin(real_window*e)/(real_window*e))
                        for w,e in zip(weights,levels) if e)
real_tail=sum(float(w) for w,e in zip(weights,levels) if e>2/real_window)
need(real_tail<=2*real_average_defect,'unitary time sinc-average tail bound')

# Divergent first moment is compatible with norm compactness/continuous time.
first_moment_controls = []
for n in (4, 16, 64, 256):
    mass = s.Rational(1, n)
    energy = n ** 2
    expectation = mass * energy
    norm_dist2 = 2 - 2*s.sqrt(1-mass)
    need(expectation == n, 'first moment diverges')
    need(norm_dist2 <= 2*mass, 'state still norm-converges to vacuum')
    first_moment_controls.append({'N': n, 'mean_energy': str(expectation),
                                  'positive_time_norm_defect_bound': str(mass)})

# Off-grade weight identity, no assumption of a homogeneous microscopic field.
Z = s.diag(1, s.I, -1, -s.I)
v = s.Matrix([s.Rational(1, 3), s.I/s.Integer(4),
              -s.Rational(1, 5), (1+s.I)/7])
grade_cases = []
for r0 in range(4):
    defect = (Z - s.I**r0*s.eye(4))*v
    lhs = s.expand((defect.H*defect)[0])
    rhs = sum(s.Abs(s.I**r-s.I**r0)**2*s.Abs(v[r])**2 for r in range(4))
    off_norm = sum(s.Abs(v[r])**2 for r in range(4) if r != r0)
    need(s.simplify(lhs-rhs) == 0, 'exact register grade identity')
    need(s.simplify(lhs-2*off_norm) >= 0, 'off-grade norm bound')
    grade_cases.append({'initial_register': r0, 'commutator_norm2': str(lhs)})

# Closed heat-filtered pair formulas and the noncommuting sharp/time limits.
r, u = s.symbols('r u', positive=True)
F = ((1+r*u)*(1-r)/((1-r*u)*(1+r)))**4
need(s.cancel(F.subs(u, 1)-1) == 0, 'finite heat time has norm1')
need(s.limit(F.subs(u, s.Rational(1, 2)), r, 1, dir='-') == 0,
     'sharp limit loses all positive-time norm')
need(s.limit(F.subs(u, s.Rational(1, 2))/(1-r)**4, r, 1, dir='-')
     == s.Rational(81, 16), 'quantitative sharp-limit rate')
need(s.cancel((-u*s.diff(F,u)).subs(u,1) + 8*r/(1-r*r)) == 0,
     'time derivative gives the energy')
heat_controls = []
for rv in (s.Rational(1, 2), s.Rational(3, 4), s.Rational(7, 8), s.Rational(15, 16)):
    exact_energy = 8*rv/(1-rv*rv)
    exact_kernel = s.cancel(F.subs({r:rv,u:s.Rational(1,2)}))
    heat_controls.append({'exp_minus_2s': str(rv), 'mean_K': str(exact_energy),
                          'F_s_log2': str(exact_kernel)})
for stime, d, t in ((0.3, math.pi, 0.4), (0.7, 0.9, 0.2)):
    n = np.arange(1, 401, dtype=float)
    q = math.exp(-2*stime)
    amplitude2 = 4*np.exp(-2*stime*n)*(1-np.cos(n*d))/n
    e_sum = float(np.sum(n*amplitude2))
    e_closed = 4*(q/(1-q)-(q*np.exp(1j*d)/(1-q*np.exp(1j*d))).real)
    b = lambda x: 4*math.log(abs(1-x*np.exp(1j*d))/(1-x))
    f_sum = math.exp(-float(np.sum(amplitude2*(1-np.exp(-t*n)))))
    f_closed = math.exp(b(q*math.exp(-t))-b(q))
    need(abs(e_sum-e_closed) < 1e-11, 'heat energy exact geometric series')
    need(abs(f_sum-f_closed) < 1e-11, 'heat time kernel exact logarithmic sum')

# Actual source energies in all four original register backgrounds.
sys.path.insert(0, str(REPO/'verification'))
src = REPO/'verification/v1033_charged_disorder.py'
spec = importlib.util.spec_from_file_location('continuous_sector_v1033', src)
mod = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = mod
spec.loader.exec_module(mod)
offset_controls = []
for nx in (8, 16, 32, 64, 128):
    energies = [mod.momentum_ground_energy(nx, k, ny=8) for k in range(4)]
    shifts = [(e-energies[1])*nx/(2*math.pi) for e in energies]
    need(abs(shifts[3]) < 1e-7, 'opposite quarter common vacuum energy')
    # Diagnostic only: the UNIFORM all-N bound is proved via BV Fourier tails.
    need(max(abs(x) for x in shifts) < 0.3, 'actual common-reference control')
    if nx == 8:
        for k in range(4):
            h = mod.qwz_cylinder(nx,8,1.,k)
            vals = np.linalg.eigvalsh(h)
            need(abs(float(np.sum(vals[vals<0]))-energies[k]) < 1e-10,
                 'Bloch and complete builder vacuum energies agree')
    offset_controls.append({'Nx':nx, 'scaled_offsets_relative_r1':shifts})
cusp = []
e0 = -np.sum(np.abs(np.linalg.eigvalsh(mod.strip_at_momentum(0.,8))))/2
for p in (1e-3, 1e-4, 1e-5):
    hp = mod.strip_at_momentum(p,8)
    ep = -np.sum(np.abs(np.linalg.eigvalsh(hp)))/2
    need(abs(np.trace(hp)) < 1e-12, 'original Bloch trace vanishes')
    cusp.append({'p':p,'one_sided_slope':float((ep-e0)/p)})
need(abs(cusp[-1]['one_sided_slope']+1) < 1e-3,
     'the two zero branches produce the stated -abs(p) cusp')

# ACTUAL nonunitary attempt: normalized exp(-s H_N) P_E Omega_N.
# Imaginary-time quadratic evolution preserves a Slater state. The Thouless
# matrix avoids cancellation of huge occupied-sea exponential factors.
actual_heat_controls = []
nx, ny = 32, 8
cols = np.arange(nx)
mask = (cols % (nx//8) >= (nx//8)//4) & (cols % (nx//8) < (nx//8)//2)
sigma = np.repeat(np.where(mask,-1.,1.), 2*ny)
for sector in (1,3):
    h = mod.qwz_cylinder(nx,ny,1.,sector)
    es, vs = np.linalg.eigh(h)
    occ, unocc = vs[:,es<0], vs[:,es>0]
    en, ep = -es[es<0]*nx/(2*math.pi), es[es>0]*nx/(2*math.pi)
    need(len(en)==len(ep)==nx*ny, 'actual quarter filling')
    aa = occ.conj().T @ (sigma[:,None]*occ)
    bb = unocc.conj().T @ (sigma[:,None]*occ)
    kk = np.linalg.solve(aa.T,bb.T).T
    loga = np.linalg.slogdet(aa)[1]
    lognorm0 = 2*loga+np.linalg.slogdet(np.eye(len(en))+kk.conj().T@kk)[1]
    need(abs(lognorm0)<1e-7, 'Thouless normalization of the actual parity state')
    step = lambda ss: np.exp(-ss*ep)[:,None]*kk*np.exp(-ss*en)[None,:]
    ss, tt = .25, .1
    ks = step(ss)
    dd = np.linalg.inv(np.eye(len(en))+ks.conj().T@ks)
    particle = ks@dd@ks.conj().T
    energy = float(np.sum(ep*particle.diagonal().real)
                   +np.sum(en*(1-dd.diagonal().real)))
    need(energy>=-1e-8,'heat-filtered excitation energy is nonnegative')
    vv = occ+unocc@ks
    cs = vv@dd@vv.conj().T
    c0 = occ@occ.conj().T
    need(np.linalg.norm(cs@cs-cs)<1e-8,'filtered full state remains a pure Slater state')
    pp = np.array([mod.sector_momentum(nx,j,sector) for j in range(nx)])
    fourier = np.exp(1j*np.outer(cols,pp))/math.sqrt(nx)
    charge_changes = []
    for residue in range(8):
        ff = fourier[:,residue::8]
        pk = np.kron(ff@ff.conj().T,np.eye(2*ny))
        dq = float(np.sum(pk.T*(cs-c0)).real)
        charge_changes.append(dq)
    need(max(abs(x) for x in charge_changes)<1e-7,
         'original heat and neutral parity preserve all eight residue charges')
    knext=step(ss+tt/2)
    ld=lambda k: np.linalg.slogdet(np.eye(len(en))+k.conj().T@k)[1]
    ftime=math.exp(ld(knext)-ld(ks))
    need(0<ftime<=1+1e-10,'actual normalized positive-time kernel')
    # Independent full Slater evolution/QR checks the Thouless chart formula.
    scaled=es*nx/(2*math.pi)
    raw=(vs.conj().T@(sigma[:,None]*occ))*np.exp(-ss*scaled)[:,None]
    qq,rr=np.linalg.qr(raw)
    cq=vs@qq@qq.conj().T@vs.conj().T
    need(np.linalg.norm(cq-cs)<1e-8,'direct full Slater covariance agrees')
    gt=qq.conj().T@(np.exp(-tt*scaled)[:,None]*qq)
    e0scaled=float(np.sum(scaled[es<0]))
    ftime_qr=math.exp(np.linalg.slogdet(gt)[1]+tt*e0scaled)
    need(abs(ftime_qr-ftime)<1e-9,'independent full Slater time kernel agrees')
    lognorm_qr=2*np.log(np.abs(rr.diagonal())).sum()+2*ss*e0scaled
    need(abs(lognorm_qr-(2*loga+ld(ks)))<1e-8,'raw sea-reference norm agrees')
    edge_changes=[]
    for y in (0,ny-1):
        i,j=2*y,2*(ny+y)
        edge_changes.append(float(np.linalg.norm((cs-c0)[i:i+2,j:j+2])))
    need(min(edge_changes)>1e-9,'actual heat filter retains responses on both edges')
    actual_heat_controls.append({'sector':sector,'Nx':nx,'Ny':ny,
        'heat_s':ss,'time_t':tt,'mean_original_H_N':energy,
        'normalized_time_kernel':ftime,'original_residue_charge_changes':charge_changes,
        'both_edge_covariance_changes':edge_changes,
        'independent_QR_kernel':ftime_qr,
        'independent_QR_covariance_error':float(np.linalg.norm(cq-cs)),
        'raw_norm_squared':math.exp(2*loga+ld(ks))})

result = {'date':'2026-10-23',
          'verdict':'SCOPED_CONTINUOUS_TIME_COMPACTNESS_EXCLUDES_NEW_SECTOR_WITH_UNCHANGED_FULL_MAP',
          'finite_checks_passed':True,
          'spectral_tail_positive_margin':str(margin),
          'real_time_average_tail_control':{'tail':real_tail,
                                            'bound':2*real_average_defect},
          'divergent_first_moment_countercontrols':first_moment_controls,
          'register_grade_controls':grade_cases,
          'heat_pair_controls':heat_controls,
          'sharp_positive_time_limit':0,'sharp_norm_at_time_zero':1,
          'source_common_reference_controls':offset_controls,
          'source_cusp_controls':cusp,
          'actual_nonunitary_QWZ_heat_filtered_pair':actual_heat_controls,
          'actual_source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),
          'general_compactness_and_uniform_reference_bounds':'written_analytic_proofs_not_extrapolations',
          'inherited_full_observable_map_required':True,
          'uniform_first_energy_moment_required':False,
          'new_P1_spinor_sector_constructed':False,
          'all_P1_no_go_claimed':False}
(OUT/'continuous-time-sector-results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'finite_checks_passed':True,'original_register_backgrounds':4,
                  'source_sizes':[8,16,32,64,128],'time_zero_norm_mismatch':True}))
