"""Transport the actual Cartan drift and test both signed native Higgs frequencies."""
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
    native = runpy.run_path(str(base/"check_tfpt_conditioned_r113_weights_20261110.py"))
    old = runpy.run_path(str(base/"check_ramond_full_return_20261110.py"))

sp,Q,m = old["sp"],old["Q"],old["m"]
M,q,V = old["M"],old["q"],old["V"]
unit = sp.eye(8)
Y6 = sp.Matrix([-2,-2,-2,3,3,0,0,0])
weak = unit[:,3]-unit[:,4]
halfswap = old["Pswap"]
color_center = sp.Matrix([1,1,-2,0,0,0,0,0])
assert halfswap*q==q
gauge_rows = [M.row(0)-M.row(1),M.row(1)-M.row(2),Y6.T*M,weak.T*M]
rows = [tuple(map(int,2*row)) for row in gauge_rows]
target = (0,0,12,-4)
q2 = tuple(map(int,2*q))
quarters = (sp.Rational(1,4),sp.Rational(3,4))


def old_omega(alpha):
    return sp.Matrix([j+alpha-4 for j in range(8)])


def signed_candidates():
    # Enumerate the possible output Ramond charges at energies 32 +/- 1/2
    # for alpha=1/4, then retain the same charge and mode at alpha=3/4.
    cap = 260
    def component_energy(j,x):
        return 8*x*x+(4*j-15)*x
    lo = [min(component_energy(j,-1),component_energy(j,1)) for j in range(8)]
    assert min(component_energy(j,x)+sum(lo)-lo[j] for j in range(8) for x in (-7,7))>cap
    options = [[(x,component_energy(j,x)) for x in (-5,-3,-1,1,3,5)] for j in range(8)]
    answer = {-1:[],1:[]}

    def rec(p,cost,w):
        j = len(p)
        if j==8:
            if (sum(p)-sum(q2))%4 or tuple(w)!=target:
                return
            r = sp.Matrix([(x-y)//2 for x,y in zip(p,q2,strict=True)])
            for sign in (-1,1):
                requested = 256+4*sign
                if sum(r)!=2*sign or cost>requested or (requested-cost)%64:
                    continue
                mode = (old_omega(quarters[0]).dot(r)-2*sign*quarters[0])/8
                assert mode.is_Integer
                assert all(old_omega(a).dot(r)-8*mode==2*sign*a for a in quarters)
                ground = (q+r).dot(q+r)/2
                assert ground+(requested-cost)//64 == 4-mode
                answer[sign].append((r,mode,ground))
            return
        for x,contrib in options[j]:
            if cost+contrib+sum(lo[j+1:])<=cap:
                rec(p+[x],cost+contrib,[z+row[j]*x for z,row in zip(w,rows,strict=True)])

    rec([],0,[0]*4)
    for sign,vals in answer.items():
        answer[sign] = sorted(vals,key=lambda v:(v[0].dot(v[0]),tuple(v[0])))
        assert len(vals)==2 and all(4-mode==ground for _,mode,ground in vals)
    return answer


def contraction(poly,charge):
    coeffs = [Q.from_sympy(x) for x in list(V.T*charge)+[q.dot(charge)]]
    value = Q.zero
    # All selected target eigenspaces are oscillator vacua.
    for word,a in poly.items():
        for _,axis in word:
            a *= -coeffs[axis]
        value += a
    return value


candidates = signed_candidates()
controls, gauge_readouts = {}, {}
time, shift = sp.Rational(1,8),sp.Rational(1,20)
gamma,mr,phi = native["gamma"],native["majorana"],native["inner"]
swap,build = native["swap"],native["native"].rep
T3 = -1j*gamma@phi(np.array([0.,1.]))
T4 = -1j*gamma@phi(np.array([0.,1j]))
H_plus = T3-1j*T4
for alpha in quarters:
    omega = old_omega(alpha)
    w = M*omega
    color = [x-sum(w[:3])/3 for x in w[:3]]
    g = (w[3]-w[4])/2
    gauge = sp.Matrix(color+[g,-g,0,0,0])
    assert g==-2 and gauge.dot(M*q)==0
    gauge_old = M.T*gauge
    remainder = w-gauge
    assert remainder[0]==remainder[1]==remainder[2] and remainder[3]==remainder[4]
    color_phases = [sp.simplify(sp.exp(4*sp.pi*sp.I*x)) for x in color]
    assert len(set(color_phases))==1 and color_phases[0]!=1
    derivative = build(0,1j*np.diag([float(g),-float(g)]),1j*np.diag(list(map(float,color))))
    G_native = -1j*(derivative+swap@derivative.conj()@swap)
    h0 = float(alpha)*gamma-1j*gamma@mr
    assert np.linalg.norm(G_native-G_native.conj().T)<1e-13
    assert np.linalg.norm(h0@G_native-G_native@h0)<1e-13
    assert np.linalg.norm(G_native@H_plus-H_plus@G_native-2*H_plus)<1e-12
    quark = np.array([x[0] in ("u","d") for x in native["native"].P]*2)
    kin = float(alpha)*gamma
    native_comm = kin@H_plus-H_plus@kin
    for sign in (-1,1):
        part = (H_plus+sign*native_comm/(2*float(alpha)))/2
        assert np.linalg.norm(part[np.ix_(quark,quark)])>1.
    gauge_readouts[str(alpha)] = {
        "SM_semisimple_drift_target":list(map(str,gauge)),
        "commuting_remainder_target":list(map(str,remainder)),
        "weak_raising_drift":str(w.dot(weak)),
        "unitary_color_cycle_at_4pi":str(color_phases[0]),
        "real_euclidean_Wilson_factor_is_retained":True,
        "native_lower_Higgs_gauge_frequency":2}
    delta = omega/8
    coeff = [Q.from_sympy(x) for x in old["cl"](V.T*delta)]
    polynomials,norms = m.code()
    images = [old["vertex_zero"](old["series_response"](f,coeff,4)) for f in polynomials]
    existing = json.loads((base/"ramond-full-return-results-20261110.json").read_text())
    scale = sp.diag(*[1/sp.sqrt(n) for n in norms])
    existing_gram = sp.Matrix(existing[str(alpha)]["full_return_Gram"])
    gram = old["cl"](scale*sp.Matrix(5,5,lambda i,j:old["inner"](images[i],images[j]))*scale)
    assert gram==existing_gram
    marked = []
    for sign in (-1,1):
        coefficients = []
        for charge,mode,ground in candidates[sign]:
            values = [contraction(x,charge) for x in images]
            coefficients.append([Q.to_sympy(x) for x in values])
            ratios = [sp.simplify(Q.to_sympy(m.conj(x)*x)/old["inner"](y,y))
                      for x,y in zip(values,images,strict=True)]
            dE = 2*sign*alpha
            gauge_charge = gauge_old.dot(charge)
            assert gauge_charge==2
            swapped = halfswap*charge
            swapped_color = color_center.dot(M*(q+swapped))%3
            assert color_center.dot(M*(q+charge))%3==0 and swapped_color in (1,2)
            assert omega.dot(swapped-charge)==16*sign
            energy = 32+dE
            ratio = ratios[0]
            returns = [ratio*sp.exp(-energy*time),ratio*sp.exp(-energy*2*time),
                       ratio*sp.exp(-energy*time-32*shift),ratio*sp.exp(-energy*(time+shift))]
            # Same fixed zero-duration unmarked q reference for all four
            # entries; extra time is assigned to the q or the marked leg.
            assert sp.exp(-(energy-gauge_charge)*(time+shift))*sp.exp(-gauge_charge*(time+shift))==sp.exp(-energy*(time+shift))
            incorrectly_unlinked_difference = ratio*sp.exp(-(energy-gauge_charge)*time)*(
                sp.exp(-32*shift)-sp.exp(-(energy-gauge_charge)*shift))
            if sign==1 and ratio:
                assert incorrectly_unlinked_difference<0<returns[2]-returns[3]
            marked.append({"sign":sign,"charge_old":list(map(str,charge)),
                           "old_mode":str(mode),"conformal_weight":str(charge.dot(charge)/2),
                           "final_unshifted_L0":str(ground),"old_intermediate_energy":str(energy),
                           "same_fixed_P0_reference_values_t_2t_shiftR_shiftH":list(map(float,returns)),
                           "P0_shift_difference":float(returns[2]-returns[3]),
                           "P0_difference_if_temporal_factor_is_omitted":float(incorrectly_unlinked_difference),
                           "old_halfswap_changes_color_center_from_0_to":int(swapped_color),
                           "old_halfswap_changes_energy_by":16*sign,
                           "exact_q_norm_ratios_for_all_five":list(map(str,ratios)),
                           "exact_ground_amplitudes_up_to_old_cocycle":[str(Q.to_sympy(x)) for x in values],
                           "rank_on_five":int(any(values))})
        assert sp.Matrix(coefficients).to_DM().convert_to(Q).rank()==2
        for row in marked[-2:]:
            row["exact_combined_same_sign_rank_on_five"] = 2
    controls[str(alpha)] = marked

print(json.dumps({
    "old_gauge_cartan_factorization_and_native_real_action":gauge_readouts,
    "actual_old_spin_and_gauge_matched_signed_frequency_controls":controls,
    "old_t":str(time),"old_positive_window_shift":str(shift),
    "physical_native_Padm_values":{"return_t":None,"return_2t":None,
                                   "ramond_shift":None,"higgs_shift":None}
},indent=2))
