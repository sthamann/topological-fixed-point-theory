"""Test the direct family lift on the original paired Pin domains."""
from contextlib import redirect_stdout
from pathlib import Path
import io
import json
import runpy
import sys

import numpy as np

repo = Path(sys.argv[1]).resolve()
base = repo/"experiments/theory-contracts/chat-triangulation-20261108"
with redirect_stdout(io.StringIO()):
    cover = runpy.run_path(str(base/"check_tfpt_cover_selector_compiler_20261110.py"))

sp, q, Q, vertex, phase = (cover[k] for k in ("sp","q","Q","vertex","phase"))
gate, signed = cover["gate"], cover["signed"]
dwell = sp.Rational(1,20)
kernel = {}
for alpha, source_rows in signed["controls"].items():
    kernel[alpha] = {}
    for sign in (-1,1):
        low,high = [row for row in source_rows if row["sign"]==sign]
        top = q+sp.Matrix(list(map(sp.sympify,low["charge_old"])))
        bottom = q+sp.Matrix(list(map(sp.sympify,high["charge_old"])))
        mid = (top+bottom)/2
        assert top-bottom==2*q and top.dot(q)==2 and mid.dot(q)==0
        for state in (top,mid,bottom):
            folded = gate["M"]@state
            assert folded[0]==folded[1]==folded[2]
            assert folded.dot(signed["Y6"])==3 and folded.dot(signed["weak"])==-1
        direct_phase = phase(-2*q,top)
        two_phase = phase(-q,top)*phase(-q,mid)
        times, partitions, two_face = {}, {}, {}
        E = sp.Rational(low["old_intermediate_energy"])
        ell = sp.Rational(1,20)
        for k in range(-3,6):
            first = vertex({():Q.one},-q,top,k)
            second = vertex(first,-q,mid,-k)
            assert not second or set(second)=={()}
            osc_coeff = Q.to_sympy(second.get((),Q.zero))
            expected = { -1:1, 0:-2, 1:1 }.get(k,0)
            assert osc_coeff==expected
            ratio = sp.simplify(two_phase*osc_coeff/(2*direct_phase))
            assert ratio=={ -1:-sp.Rational(1,2),0:1,1:-sp.Rational(1,2) }.get(k,0)
            times[k] = sp.simplify(ratio*sp.exp(-8*k*dwell))
            x = sp.exp(-ell*E)
            xp = x*sp.exp(-8*k*ell)
            first_graph_factor = sp.sqrt((1+xp*xp)/(1+x*x))
            second_graph_factor = sp.sqrt((1+x*x)/(1+xp*xp))
            assert sp.simplify(first_graph_factor*second_graph_factor)==1
            if k in (-1,0,1):
                two_face[str(k)] = [str(sp.exp(-8*k*ell)),str(sp.exp(8*k*ell))]
            if k in (1,2):
                contributions = {}
                for word,value in first.items():
                    part = vertex({word:value},-q,mid,-k).get((),Q.zero)
                    c = sp.simplify(two_phase*Q.to_sympy(part)/(2*direct_phase))
                    contributions["+".join(str(n) for n,_ in word)] = str(c)
                assert sum(map(sp.Rational,contributions.values()))==ratio
                partitions[str(k)] = contributions
        assert partitions["1"]=={"2":"1/2","1+1":"-1"}
        assert set(partitions["2"].values())=={"1/3","-1","2/3"}
        summed = 1-sp.cosh(8*dwell)
        assert sp.simplify(sum(times.values())-summed)==0
        direct = cover["old"][alpha][str(sign)]["direct_old_mark"]
        kernel[alpha][str(sign)] = {
            "common_fixed_old_time": "1/8",
            "intermediate_dwell": str(dwell),
            "whole_vertex_same_color_and_reference":direct,
            "split_relative_factors_for_first_mode_k":{str(k):str(times[k]) for k in (-1,0,1)},
            "second_face_relative_mode_transport_on_the_two_halves":two_face,
            "old_normalized_graph_factors_cancel_in_complete_word":True,
            "individual_old_oscillator_paths_before_time_factor":partitions,
            "split_absolute_complex_pairs":{str(k):[float(times[k])*v for v in direct] for k in (-1,0,1)},
            "split_sum_relative_to_whole":str(summed),
            "split_sum_relative_float":float(summed),
            "split_with_old_R_inner_filter": [0.,0.],
        }

# In the original oscillator expansion the first step creates level k+1.
# Full contraction into the final vacuum has generating function (1-u)^2.
u = sp.Symbol("u")
analytic_contraction = sp.exp(-2*sum(u**j/sp.Integer(j) for j in range(1,8))).series(u,0,8).removeO()
assert sp.expand(analytic_contraction-(1-u)**2)==0
d = sp.Symbol("d",real=True)
assert sp.limit((1-sp.cosh(8*d))/d**2,d,0)==-32

native = gate["native"]
rep, gamma, phi = (native[k] for k in ("native","gamma","inner"))
mass = native["majorana"]
marks = [phi(np.array(v)) for v in ((1,0),(1j,0),(0,1),(0,1j))]
tau1,tau2,tau3 = (cover[k] for k in ("sigma1","sigma2","sigma3"))
Id = np.eye(96)
P = [np.array(x,dtype=complex) for x in (tau1,tau2,tau3)]
pin_su2 = [np.kron(P[0],gamma),np.kron(P[1],Id),np.kron(P[2],gamma)]
family_probes = [np.kron(P[1],v) for v in marks+[mass]]
for S in pin_su2:
    assert np.array_equal(S@S,np.eye(192))
    assert all(np.linalg.norm(S@v-v@S)<1e-13 for v in family_probes)
assert np.array_equal(pin_su2[0]@pin_su2[1]-pin_su2[1]@pin_su2[0],2j*pin_su2[2])
N = (pin_su2[0]+1j*pin_su2[1])/2
assert np.linalg.norm(N)>1 and np.linalg.norm(N@N)<1e-13
# The simple swap is a spin-half action, but the original paired free
# CAR admits a spin-one two-mode sector. Test it using independent
# positive-energy electron and conjugate modes, not 96 new particles.
C = np.kron(np.eye(2),native["swap"])
assert all(np.linalg.norm(C@S.conj()@C+S)<1e-13 for S in pin_su2)
phi_background = mass+.55*marks[2]
K_zero = -np.kron(P[1],phi_background)
assert np.linalg.norm(C@K_zero.conj()@C+K_zero)<1e-13
electron = [i for i,row in enumerate(rep.P) if row[0]=="e"]
labels = [electron,[i+rep.NP for i in electron]]
zero2,zero3 = np.zeros((2,2)),np.zeros((3,3))


def real_gauge(a):
    j = native["swap"]
    return -1j*(a+j@a.conj()@j)


em = real_gauge(rep.rep(0,1j*np.diag([1,-1])/2,zero3))
em += real_gauge(rep.rep(3j,zero2,-1j*np.eye(3)))/6
for charge,ids in zip((-1,1),labels,strict=True):
    np.testing.assert_allclose(em[np.ix_(ids,ids)],charge*np.eye(len(ids)),atol=1e-14)
assert np.linalg.norm(em@phi_background-phi_background@em)<1e-13
copy_energies = []
for ids in labels:
    selected = np.array(ids+[i+96 for i in ids])
    restricted = K_zero[np.ix_(selected,selected)]
    energies, U = np.linalg.eigh(restricted)
    largest = energies[-1]
    eigspace = U[:,abs(energies-largest)<1e-10]
    assert eigspace.shape[1]==2 and largest>0
    spin3 = pin_su2[2][np.ix_(selected,selected)]
    eigs, vectors = np.linalg.eigh(eigspace.conj().T@spin3@eigspace)
    np.testing.assert_allclose(eigs,[-1,1],atol=1e-13)
    top = eigspace@vectors[:,1]
    minus = N.conj().T[np.ix_(selected,selected)]@top
    assert abs(np.linalg.norm(minus)-1)<1e-13
    assert np.linalg.norm(restricted@minus-largest*minus)<1e-13
    copy_energies.append(float(largest))
np.testing.assert_allclose(copy_energies[0],copy_energies[1],atol=1e-13)
lower2 = sp.Matrix([[0,0],[1,0]])
Jm = sp.kronecker_product(lower2,sp.eye(2))+sp.kronecker_product(sp.eye(2),lower2)
top2 = sp.eye(4)[:,0]
middle2 = Jm*top2/sp.sqrt(2)
bottom2 = Jm*middle2/sp.sqrt(2)
singlet = (sp.eye(4)[:,1]-sp.eye(4)[:,2])/sp.sqrt(2)
assert Jm*bottom2==Jm*singlet==sp.zeros(4,1)
assert (Jm**2/2)*top2==bottom2
Jz2 = sp.diag(1,0,0,-1)
assert [int((v.H*(sp.I*sp.pi*Jz2).exp()*v)[0]) for v in (top2,middle2,bottom2)]==[-1,1,-1]


# Fixed SM and grading give one generation matrix per isotypic block.
# Use both actual Yukawa squares on each left doublet to remove every
# non-scalar generation matrix. Then invertibility forces the right side.
def exact_matrix(a):
    return sp.Matrix([[sp.Rational(float(x.real))+sp.I*sp.Rational(float(x.imag)) for x in row] for row in a])


y = {}
for species in ("u","d","nu","e"):
    left,right = [[rep.IDX[(species,side,0,i)] for i in range(3)] for side in ("L","R")]
    y[species] = exact_matrix(marks[2][np.ix_(left,right)])
    assert y[species].det()!=0
generation_commutant = {}
for label,f,g in (("quark","u","d"),("lepton","nu","e")):
    I = sp.eye(3)
    first, second = ((y[k]*y[k].H).applyfunc(sp.expand) for k in (f,g))
    matrices = [sp.kronecker_product(a,I)-sp.kronecker_product(I,a.T) for a in (first,second)]
    system = matrices[0].col_join(matrices[1])
    rank = system.to_DM().convert_to(sp.QQ_I).rank()
    assert rank==8
    generation_commutant[label]=9-rank
electron_masses = np.linalg.svd(.55*marks[2][np.ix_(
    [rep.IDX[("e","L",0,k)] for k in range(3)],
    [rep.IDX[("e","R",0,k)] for k in range(3)])],compute_uv=False)
assert min(abs(a-b) for i,a in enumerate(electron_masses) for b in electron_masses[:i])>1e-8
nu_R = [rep.IDX[("nu","R",0,k)] for k in range(3)]
assert exact_matrix(mass[np.ix_(nu_R,[i+rep.NP for i in nu_R])]).det()!=0

# Exact frequency check for the unchanged flat massive normal generator.
# Its energies are +/-sqrt((n*pi/2)^2+lambda_a), with algebraic lambda_a.
# Nonzero algebraic differences between two energies with at least one
# nonzero n would give a nonzero polynomial for pi^2 over algebraic numbers.
# At n=m=0 check the only possible shift 8 exactly on the mass polynomials.
assert gate["bound"]<64 and 4*gate["bound"]<16**2
z = gate["var"]
F = sp.Poly.from_list([1],z,domain=sp.QQ)
for poly in gate["polynomials"]:
    assert poly.eval(0)!=0
    coefficients = poly.all_coeffs()
    even = [v for c in coefficients for v in (c,0)][:-1]
    F *= sp.Poly.from_list(even,z,domain=sp.QQ)
assert F.degree()==30
prime_certificate = None
for p in (65537,65539,65543):
    reduced = gate["modular"](F,p)
    if reduced.gcd(reduced.shift(8)).degree()==0:
        prime_certificate=p
        break
assert prime_certificate is not None

# For any real constant omega, the only normal-local multiplication
# candidate with [K,X]=omega X and the given wall conditions has
# A(w)=A0 cos(omega(w+1)), B(w)=-i A0 sin(omega(w+1)).
# Check the two coefficient equations before using the full-rank Phi.
w, omega = sp.symbols("w omega",real=True)
A = sp.cos(omega*(w+1))
B = -sp.I*sp.sin(omega*(w+1))
assert sp.simplify(sp.I*B.diff(w)-omega*A)==0
assert sp.simplify(sp.I*A.diff(w)-omega*B)==0
assert B.subs(w,-1)==0

print(json.dumps({
    "actual_old_whole_and_split_local_words":kernel,
    "whole_vs_split_contact":{
        "contraction_series":"(1-u)^2",
        "only_contributing_mode_pairs_first_k":[-1,0,1],
        "zero_mode_split_equals_whole":True,
        "coincident_unsmeared_product_over_whole":0,
        "leading_old_cylinder_contact_coefficient_for_L_minus_k_L_k_over_2": "-1/2",
    },
    "native_paired_Pin_family_zero_mode":{
        "internal_common_generation_commutants":generation_commutant,
        "gauge_grading_and_field_compatible_internal_scalar_blocks":["quark","antiquark","lepton_plus_antilepton"],
        "full_normal_multiplication_commutant":"M2(C) direct_sum M2(C) direct_sum M2(C)",
        "has_time_compatible_pin_spin_half":True,
        "contains_original_spin_one_ladder":False,
    },
    "native_conditional_free_CAR_two_mode_family":{
        "number_of_actual_positive_energy_modes":4,
        "common_one_particle_energy_per_doublet":copy_energies,
        "electron_and_antiparticle_combined_electromagnetic_charge":0,
        "two_particle_sector":"spin 1 plus singlet",
        "max_spin_one_copies_at_one_energy_in_entire_free_electron_positron_pair_channel":3,
        "three_step_grade":[-1,1,-1],
        "same_fixed_pin_loop_phase_squared_is_scalar":-1,
        "actual_old_Ramond_endpoint_identification":None,
    },
    "native_paired_Pin_local_modes":{
        "nonzero_normal_local_single_frequency_multiplier_with_full_rank_Phi":False,
        "no_energy_difference_8k_for_nonzero_integer_k_on_unchanged_flat_chart":True,
        "exact_mass_shift_modular_prime":prime_certificate,
        "no_infinite_arithmetic_ladder_at_any_fixed_time_calibration_on_flat_massive_one_particle_chart":True,
    },
    "physical_full_Padm_same_source":None,
},indent=2))
