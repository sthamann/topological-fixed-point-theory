"""Common-CAR state transport and local/reflected source tests on the actual Pin pair."""

import contextlib
from fractions import Fraction
import io
import json
from pathlib import Path
import runpy
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad


repo = Path(sys.argv[1]).resolve()
base = repo/"experiments/theory-contracts/chat-triangulation-20261108"
with contextlib.redirect_stdout(io.StringIO()):
    native = runpy.run_path(str(base/"check_tfpt_conditioned_r113_weights_20261110.py"))
m = native["majorana"]
v = native["inner"](np.array([0.,1.]))
phi = m+.55*v
swap = native["swap"]
eye = np.eye(96)
z = np.zeros_like(m)
sheet_sign = np.block([[eye,z],[z,-eye]])
j = np.block([[swap,z],[z,swap]])
for op in (m,v):
    np.testing.assert_allclose(swap@op.conj()@swap,op,atol=1e-14)


def operator(k, field):
    return np.block([[k*eye,field],[field,-k*eye]])


def rotation(k, field):
    values,u = np.linalg.eigh(field)
    half = .5*np.arctan(values/k)
    c = (u*np.cos(half))@u.conj().T
    s = (u*np.sin(half))@u.conj().T
    return np.block([[c,-s],[s,c]])


def projector(k,field):
    values,u = np.linalg.eigh(operator(k,field))
    occupied = u[:,values<0]
    return occupied@occupied.conj().T


transports = {}
for k in (np.pi/2,3*np.pi/2,-np.pi/2,-3*np.pi/2):
    u = rotation(k,phi)@rotation(k,m).conj().T
    inverse_sector = rotation(-k,phi)@rotation(-k,m).conj().T
    norm_error = np.max(abs(u@u.conj().T-np.eye(len(u))))
    state_error = np.max(abs(u@projector(k,m)@u.conj().T-projector(k,phi)))
    reflection_error = np.max(abs(sheet_sign@u@sheet_sign-inverse_sector))
    real_error = np.max(abs(j@u.conj()@j-u))
    time_error = np.linalg.norm(operator(k,phi)@u-u@operator(k,m))
    assert max(norm_error,state_error,reflection_error,real_error)<1e-12
    assert time_error>1.
    transports[str(k)] = {"CAR_unitarity_error":float(norm_error),
                          "state_intertwiner_error":float(state_error),
                          "Green_reflection_covariance_error":float(reflection_error),
                          "D_sea_real_reflection_error":float(real_error),
                          "fixed_time_intertwiner_HS_defect":float(time_error)}

# Choose an actual native quark mass eigenmode; M_R vanishes there.
quark = np.array([f[0] in ("u","d") for f in native["native"].P]*2)
eigenmass = float(np.linalg.eigvalsh(phi[np.ix_(quark,quark)])[-1])
nodes,weights = leggauss(120)


def smeared(k,a,b):
    x = (a+b)/2+(b-a)*nodes/2
    dx = (b-a)*weights/2
    bump = np.exp(1-1/(1-nodes**2))
    return np.exp(-1j*k[:,None]*x)@(dx*bump)/np.sqrt(np.sum(dx*bump*bump))


locality = {}
for cutoff in (32,64,128):
    k = 2*np.pi*(np.arange(-cutoff,cutoff+1)+.25)
    positive = smeared(k,.12,.29)
    negative = smeared(k,-.32,-.17)
    result = np.vdot(negative,np.sin(.5*np.arctan(eigenmass/k))*positive)
    locality[str(cutoff)] = [result.real,result.imag]
np.testing.assert_allclose(locality["64"],locality["128"],atol=2e-8)
assert np.linalg.norm(locality["128"])>.01

# The old q pure states are not quasifree. Charge fixes the total occupation
# on the active band exactly; the Wick state from their covariance would
# have independent Bernoulli variance sum p(1-p).
packet_cov = json.loads((base/"tfpt-ramond-covariance-and-cover-results-20261110.json").read_text())
wick = {}
for quarter,case in packet_cov["cases"].items():
    values = []
    for packet in case["q_output"]:
        diagonal = list(map(Fraction,packet["N_diagonal_by_channel_then_mode"]))
        variance = sum(p*(1-p) for p in diagonal)
        assert variance>0 and packet["F_identically_zero"]
        values.append({"exact_actual_variance":0,"Wick_variance_exact":str(variance),
                       "Wick_variance":float(variance),
                       "fractional_modes":sum(0<x<1 for x in diagonal)})
    wick[quarter]=values

# Normal zero branch of the actual flat paired K. When V marks both sides,
# its even sheet matrix element averages the +/- native mass energies.
eigs,u = np.linalg.eigh(phi)
higgs_weights = np.sum(abs(u.conj().T@v)**2,axis=1)/len(v)


def kernel(u_time, mass):
    def integrand(mu):
        e = np.hypot(mass,mu)
        return (np.exp(-u_time*e)+np.exp(-(2-u_time)*e))/(1+np.exp(-2*e))
    return quad(integrand,0,np.inf,epsabs=2e-11)[0]/(2*np.pi)


def halfword(u_time):
    return sum(w*kernel(u_time,mass) for w,mass in zip(higgs_weights,eigs,strict=True))


def direct_fourier(u_time):
    k = np.pi*(np.arange(32768)+.5)
    mass2 = eigs**2
    e = np.hypot(k[:,None],eigs)
    difference = -mass2/(e*(e+k[:,None]))
    weights = np.sin(k*u_time)/2
    leading = higgs_weights.sum()/(4*np.sin(np.pi*u_time/2))
    # Discrete summation by parts bounds the uncomputed oscillatory tail.
    result = leading+weights@(difference@higgs_weights)
    tail_bound = np.sum(higgs_weights*mass2)/(4*k[-1]**2*abs(np.sin(np.pi*u_time/2)))
    return result,tail_bound


times = (.1,.3,.6)
gram = np.array([[halfword(a+b) for b in times] for a in times])
assert np.linalg.eigvalsh(gram)[0]>0
raw = {}
for u_time in (.3,.5,.7):
    a = halfword(u_time)
    f,bound = direct_fourier(u_time)
    assert abs(a-f)<bound+1e-10
    raw[str(u_time)] = {"spectral_integral":a,"Fourier_reference":float(f),"tail_bound":float(bound)}
np.testing.assert_allclose(kernel(.5,0.),1/(4*np.sin(np.pi/4)),atol=1e-13)

# Even neutral test under the unmodified old Green reflection. The two
# native neutral Majorana eigenchannels are orthogonal; use real temporal
# test functions. The diagonal spatial sheet sees the sum of both energies.
positive_times = .205+.085*nodes
dx = .085*weights
bump = np.exp(1-1/(1-nodes**2))
measure = dx*bump/np.sqrt(np.sum(dx*bump*bump))


def smeared_cross(mass):
    def integrand(mu):
        e = np.hypot(mass,mu)
        return ((measure@np.exp(-e*positive_times))**2
                +(measure@np.exp(-e*(1-positive_times)))**2)/(1+np.exp(-2*e))
    return quad(integrand,0,np.inf,epsabs=2e-13)[0]/(2*np.pi)


neutral_masses = np.linalg.eigvalsh(m)
neutral_masses = neutral_masses[neutral_masses>1e-9][:2]
neutral_legs = [smeared_cross(mass) for mass in neutral_masses]
neutral_expected = -neutral_legs[0]*neutral_legs[1]
# Independent finite CAR matrix expectation on the exact 4-mode
# compression C. Modes are (+,-) on each of the two channels.
c = np.eye(4,dtype=complex)/2
for index,x in enumerate(neutral_legs):
    c[2*index+1,2*index] = 1j*x
    c[2*index,2*index+1] = -1j*x
p,u = np.linalg.eigh(c)
annihilation = np.zeros((4,16,16),complex)
for i in range(4):
    for bits in range(16):
        if bits>>i&1:
            annihilation[i,bits^(1<<i),bits] = (-1)**((bits&((1<<i)-1)).bit_count())
field = np.einsum("ij,jab->iab",u,annihilation)
word = field[3].conj().T@field[1]@field[0].conj().T@field[2]
diagonal_weight = [np.prod([p[i] if bits>>i&1 else 1-p[i] for i in range(4)])
                   for bits in range(16)]
neutral_direct = np.dot(diagonal_weight,word.diagonal())
np.testing.assert_allclose(neutral_direct,neutral_expected,atol=1e-13)
assert neutral_direct.real<0
num_plus = field[0].conj().T@field[0]-np.eye(16)/2
num_minus = field[1].conj().T@field[1]-np.eye(16)/2
centered_density = np.dot(diagonal_weight,(num_minus@num_plus).diagonal())
np.testing.assert_allclose(centered_density,-neutral_legs[0]**2,atol=1e-13)

# Evaluate the same two-channel bilinear with the existing ordered NS Green
# functional on the four-label cover. The normal zero mode has spatial
# generator -b*tau2*Phi. A single sheet sees both energies with weight 1/2.
def ordered_green_leg(mass):
    forward = measure@np.exp(-mass*positive_times)
    around_cover = measure@np.exp(-mass*(1-positive_times))
    return (forward**2+around_cover**2)/(2*(1+np.exp(-2*mass)))


green_legs = [ordered_green_leg(mass) for mass in neutral_masses]
green_return = green_legs[0]*green_legs[1]
assert green_return>0

print(json.dumps({"native_CAR_transport":transports,
                  "cross_reflection_support_native_quark_mass":eigenmass,
                  "disjoint_half_support_matrix_element":locality,
                  "q_packets_Wick_failure":wick,
                  "phase_adjusted_pure_D_sea_normal_zero_Higgs_kernel":raw,
                  "positive_two_leg_kernel_three_time_gram_min":float(np.linalg.eigvalsh(gram)[0]),
                  "positive_two_leg_kernel_three_time_gram_rank":int(np.linalg.matrix_rank(gram)),
                  "unmodified_Green_reflection_even_neutral_word":{
                      "native_neutral_masses":neutral_masses.tolist(),
                      "positive_phase_removed_legs":neutral_legs,
                      "Wick_return":float(neutral_expected),
                      "independent_finite_CAR_return":float(neutral_direct.real),
                      "centered_density_return":float(centered_density.real),
                      "reflection_convention":"old charge exchange with minus on each field and reversed word"},
                  "same_bilinear_ordered_free_NS_Green_comparison":{
                      "charge_reflected_legs":[float(x) for x in green_legs],
                      "reflected_return":float(green_return),
                      "scope":"Same masses, test supports and reflection on the four-label cover. Free ordered Green functional, not the interacting P1 source."},
                  "scope":"Actual flat static paired Pin operators and old q packets. Phase-adjusted two-leg form only. No mixed physical P1 histories, actual QTheta or R113 producer."},
                 indent=2))
