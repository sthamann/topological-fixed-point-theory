"""Original F1 source-first gate: R113 empty cap, actual P, Higgs kinetic
Dirichlet form, hidden generator moments and nontracial modular clock.
This is a necessary SOURCE test, not an artificial Q0 packet or F1 proof.
"""
from pathlib import Path
import json
import sympy as sp

ROOT=Path(__file__).resolve().parents[2]
NATIVE=ROOT/"experiments/theory-contracts/original-paper-source-reconciliation-20261102/contribution_20261106Q/domain/native-diagnostic-v1.json"
R113=ROOT/"experiments/theory-contracts/chat-triangulation-20261108/R113-ORIGINAL-UND-REFERENZ-20261109.txt"
KINETIC=ROOT/"experiments/theory-contracts/chat-triangulation-20261108/tfpt-common-kinetic-clock-rule-results-20261111.json"
METROLOGY=ROOT/"documents/archive--tfpt-45--source_extracts--05_metrology_source/sections/002-subsection-closed-texorpdfstring-41-41-chain-after-famil.tex"
SPATIAL=ROOT/"research/originlab/F1-NATIVE-SPATIAL-RETURN-20261115.md"
PHASESCAN=ROOT/"research/originlab/F1-R113-PHASE-SCAN-20261115.md"


def require(ok,msg):
    if not ok:
        raise ArithmeticError(msg)


def original_packet():
    native=json.loads(NATIVE.read_text(encoding="utf-8"))
    r113=R113.read_text(encoding="utf-8")
    scan=PHASESCAN.read_text(encoding="utf-8")
    # A whole source bundle may appear later; report it for review instead of failing.
    bundle_present=native.get("whole_source_bundle") is not None
    require(native["coordinate_role"]=="spatial_tangential_circle_NOT_OS_time",
            "original source time changed")
    require("C_{j,t}^{00}(b_out,b_in)" in r113 and
            "<Ω_F,Γ(R_j[b])Ω_F>=1" in r113 and
            "P_adm Q_A[B] P_adm" in r113,
            "original vacuum and outer physical selector contract changed")
    require("centre projected trace" in " ".join(scan.split()) and
            "whole-source R113 operation" in scan,
            "conditional phase-scan source distinction changed")
    return {
      "domain":native["domain_id"],
      "full_original_source_packet_present":bundle_present,
      "source_packet_review_required":bundle_present,
      "original_spatial_coordinate_is_OS_time":False,
      "original_source_operator":"Q[b]=s[b]*Gamma(R[b])",
      "original_pre_projector_zero_fermion_cap":"1",
      "original_empty_boson_endpoint":"I(w[b]*s[b])",
      "outer_physical_selector_needed":True,
      "conditional_center_trace_is_original_empty_cap":False,
    }


def original_empty_cap_projector():
    r,s=sp.symbols("r s")
    T=sp.diag(1,r,s,r*s)  # Fock order: empty, first, second, pair
    vac=sp.eye(4)[:,0]
    pair=sp.eye(4)[:,3]
    parity=sp.diag(1,-1,-1,1)
    number=sp.diag(0,1,1,2)
    v=(vac+pair)/sp.sqrt(2)
    P=v*v.H
    c=sp.simplify((vac.H*P*T*P*vac)[0])
    require((vac.H*T*vac)[0]==1 and
            sp.simplify(sp.trace(T)-(1+r)*(1+s))==0 and
            P*P==P and P.H==P and P*parity==parity*P and
            P*number!=number*P and c==(1+r*s)/4,
            "original zero-Fock cap and even physical selector conflated")
    P0=sp.diag(1,0,0,0)
    require(P0*number==number*P0 and (vac.H*P0*T*P0*vac)[0]==1,
            "number preserving vacuum selector changed")
    z=sp.symbols("z",nonzero=True)
    Tz=sp.diag(1,r,s*z,r*s*z)
    make_pair=pair*vac.H
    break_pair=make_pair.H
    unmarked=(vac.H*Tz*vac)[0]
    marked=(vac.H*break_pair*Tz*make_pair*vac)[0]
    require(unmarked==1 and marked==r*s*z and
            sp.diff(unmarked,z)==0 and sp.diff(marked,z)==r*s,
            "the original empty cap was silently replaced by a trace")
    return {
      "empty_original_fermion_cap":"1 for any R",
      "full_Fock_trace_not_same_source":"(1+r)*(1+s)",
      "allowed_even_parity_projected_cap":"(1+r*s)/4",
      "parity_preserved_but_number_mixed":True,
      "number_preserving_vacuum_retaining_cap":"1",
      "even_ordered_pair_source_independent_from_flat_unmarked_phase_scan":True,
      "unmarked_empty_cap_phase_scan":"1",
      "same_Fock_source_ordered_even_pair":"r*s*z",
      "all_original_R113_marks_follow_from_unmarked_phasescan_alone":False,
      "physical_original_P_adm_number_action_computed":False,
    }


def original_higgs_kinetic_form():
    data=json.loads(KINETIC.read_text(encoding="utf-8"))
    paper=METROLOGY.read_text(encoding="utf-8")
    require("Z_H=1/(16\\pi^2)" in paper and
            data["U7_on_original_profile"]=="13*a/100",
            "original kinetic normalization or U7 profile changed")
    all_values={}
    for sign in ("-1","1"):
        row=data["kinetic_coefficients_in_units_one_over_pi_R_Z"][sign]
        square=sum(sp.Rational(g[2])**2/(4*sp.Integer(g[3]))
                   for g in row["gradient_components"])
        require(len(row["gradient_components"])==9 and
                square==sp.Rational(5117,5000) and
                square==sp.Rational(row["whole_local_gradient"]) and
                square/sp.Rational(row["one_moving_mode"])==
                    sp.Rational(40936,169),
                "original entire-field gradient did not reproduce 9 probes")
        all_values[sign]=str(square)
    R=sp.symbols("R",positive=True)
    Z=1/(16*sp.pi**2)
    full=sp.simplify(sp.Rational(5117,5000)/(sp.pi*R*Z))
    require(sp.simplify(full-sp.Rational(10234,625)*sp.pi/R)==0,
            "original local kinetic scale changed")
    return {
       "printed_original_Higgs_residue":"Z_H=1/(16*pi^2)",
       "original_U7_readout":"13*a/100",
       "nine_original_nonzero_gradient_directions":9,
       "source_gradient_coefficients_full":all_values,
       "full_source_over_one_mode":"40936/169",
       "full_source_kinetic_local_coefficient":"(10234/625)*pi/R",
       "physical_U7_vacuum_expectation_derived":False,
       "physical_source_P_adm_reduced_energy_form_derived":False,
       "physical_OS_time_generator_derived_from_local_Z_H":False,
    }


def exact_generator_moment_defect():
    # For true positive-time H and observed source word vectors B:
    # S=M2-M1 M0^-1 M1=B*H(I-P_B)H B >=0.
    H=sp.diag(0,1,2)
    B=sp.Matrix([[1,0],[0,1],[1,0]])
    m0=B.H*B
    m1=B.H*H*B
    m2=B.H*H*H*B
    proj=B*m0.inv()*B.H
    defect=m2-m1*m0.inv()*m1
    require(defect==B.H*H*(sp.eye(3)-proj)*H*B
            and defect==sp.diag(2,0) and defect.rank()==1,
            "moment hidden-source rank/PSD identity failed")
    spatial=SPATIAL.read_text(encoding="utf-8")
    require("192(k+1)" in spatial and
            "invertible products" in spatial and
            "spatial" in spatial.lower(),
            "native full96 infinite spatial memory evidence changed")
    return {
      "common_physical_generator_closure_test":
        "M2-M1*M0^-1*M1=B*H(I-P_B)H B >=0",
      "finite_exact_comparison_hidden_rank":1,
      "native_spatial_full96_first_hidden_rank":192,
      "native_spatial_Hankel_rank":"192*(k+1)",
      "spatial_full96_return_proves_physical_OS_time":False,
      "original_common_Higgs_family_clock_stress_M0_M1_M2_available":False,
    }


def tracial_modular_time_obstruction():
    e01=sp.eye(3)[:,0]*sp.eye(3)[1,:]
    def deriv(rho):
        lr=sp.diag(*[sp.log(rho[k,k]) for k in range(3)])
        return lr*e01-e01*lr
    source_tracial=sp.eye(3)/3
    nontracial=sp.diag(sp.Rational(1,2),sp.Rational(1,3),
                       sp.Rational(1,6))
    require(deriv(source_tracial)==sp.zeros(3) and
            sp.simplify(deriv(nontracial)[0,1]-
                        sp.log(sp.Rational(3,2)))==0,
            "finite source modular flow state dependence changed")
    return {
      "M3_symmetric_source_state":"I3/3",
      "tracial_state_modular_Hamiltonian":"scalar; modular automorphism trivial",
      "tracial_M3_modular_flow_nontrivial":False,
      "positive_nontracial_example":"diag(1/2,1/3,1/6)",
      "nontrivial_modular_E01_frequency":"log(3/2)",
      "physical_source_nontracial_state_derived":False,
    }


def controls():
    return {
      "original_native_status":original_packet(),
      "same_source_empty_cap_and_selector":original_empty_cap_projector(),
      "native_full_Higgs_kinetic_direction":original_higgs_kinetic_form(),
      "hidden_source_generator_moments":exact_generator_moment_defect(),
      "state_and_modular_clock":tracial_modular_time_obstruction(),
      "original_full_R113_four_numbers_computed":False,
      "original_common_positive_generator_constructed":False,
      "original_physical_3plus1D_E8_SM_spin2_derived":False,
      "F1_solved":False,
    }


if __name__=="__main__":
    print(json.dumps(controls(),indent=2,sort_keys=True))
