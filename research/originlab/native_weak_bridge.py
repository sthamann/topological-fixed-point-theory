"""Test the split weak word on the existing full-96 R113 field chart."""

import json
from collections import Counter
from pathlib import Path
import subprocess
import sys

from coherence import R104_SOURCES, car_word, det_word, inner, r104, thermal
from source import MARKS, NATIVE
from verify import digest
from weak_source import native_d8_legs

RUNTIME = NATIVE.parents[1]/"runtime"
COMPARATOR = Path("experiments/theory-contracts/chat-triangulation-20261108")
INPUTS = (RUNTIME/"assemble.py", RUNTIME/"repo_paths.py", COMPARATOR/"check_mixed_higgs_eich.py",
          Path("verification/v252_full_finite_triple.py"), *R104_SOURCES)


def analyze(domain):
    import numpy as np
    from assemble import Assembler, HIGGS
    from check_mixed_higgs_eich import gauge_operators

    assembler = Assembler(domain)
    model = assembler.original
    weak, hypercharge = gauge_operators()
    n = model.NH
    zero2, zero3 = np.zeros((2,2),complex), np.zeros((3,3),complex)
    swap = model.Kswap()

    def gauge_generator(q, m):
        t = model.rep(0,1j*q,1j*m)
        return -1j*(t+swap@t.conj()@swap)

    colour_cartan = [gauge_generator(zero2,np.diag(h)) for h in ((-1,1,0),(0,-1,1))]
    weights = [6*hypercharge,2*weak[2],*colour_cartan]
    integer_weights = np.array([np.rint(h.diagonal().real).astype(int) for h in weights]).T
    for index, (h, q) in enumerate(zip(weights, integer_weights.T)):
        # 6Y may require one rounding because Y is stored in binary64.
        tolerance = 1e-14 if index==0 else 0
        if not np.allclose(h,np.diag(q),atol=tolerance,rtol=0):
            raise ArithmeticError("unexpected exact gauge charge dictionary")

    def projector(weight):
        return np.diag(np.all(integer_weights==weight,axis=1).astype(int))

    p0 = projector((0,0,0,0))
    anti_sterile = [model.NP+model.IDX[("nu","R",0,g)] for g in range(model.NG)]
    p_sterile = np.diag(np.isin(np.arange(n),anti_sterile).astype(int))
    lepton = [[model.NP+model.IDX[(flavour,"L",0,g)] for g in range(model.NG)]
              for flavour in ("e","nu")]
    down = [[model.IDX[("d","R",colour,g)] for g in range(model.NG)] for colour in range(3)]
    quark = [[[model.IDX[(flavour,"L",colour,g)] for g in range(model.NG)] for flavour in ("u","d")]
             for colour in range(3)]
    grade = integer_weights[:,0] % 3
    cross = grade[:,None]!=grade[None,:]
    generators = [assembler.gamma,assembler.DM,*assembler.V,*weak,hypercharge,*colour_cartan]
    if any(np.count_nonzero(letter[cross]) for letter in generators):
        raise ArithmeticError("native generators do not preserve all colour centre blocks")
    # Rebuild each 3x3 block independently over the source's decimal
    # rationals, then check the generated matrices against those blocks.
    import sympy as sp

    def coupling(name):
        return sp.Matrix([[sp.Rational(str(re))+sp.I*sp.Rational(str(im)) for re,im in row]
                          for row in assembler.domain["couplings"][name]])

    def nonzero_block(matrix, rows, cols, expected):
        actual = matrix[np.ix_(rows,cols)]
        if not np.allclose(actual, np.asarray(expected,dtype=complex),rtol=0,atol=1e-14):
            raise ArithmeticError("recovered Higgs block does not match its original family matrix")
        return expected.rank()

    yukawas = {name:coupling(name) for name in ("Yu","Yd","Ye","Ynu")}
    ynu, yd = yukawas["Ynu"],yukawas["Yd"]
    ranks = []
    string_checks = []
    polar_checks = []
    components = native_d8_legs()
    for i in range(3):
        colour_ranks = []
        for j in range(2):
            left,middle,right = components[2*i+j]
            minus_left = tuple(-k for k in left)
            early_weight, late_weight, middle_weight = map(_root_weight,(right,minus_left,middle))
            t3 = early_weight[1]
            p_in = projector(early_weight)
            p_out = projector(late_weight)
            p_gamma = projector(middle_weight)
            # The unmodified Higgs matrices are derivative directions of
            # a real quaternion, so resolve their conjugate complex blocks.
            first = 2*j
            vertex = (assembler.V[first]-1j*assembler.V[first+1])/2
            if np.count_nonzero(vertex*(integer_weights[:,0,None]-integer_weights[None,:,0]!=3)) or np.count_nonzero(
                    vertex*(integer_weights[:,1,None]-integer_weights[None,:,1]!=t3)):
                raise ArithmeticError("the extracted block is not the original charged Higgs component")
            a = p_in@vertex@p0
            delta = p_out@vertex@p_gamma
            sign = -1 if j==0 else 1
            row = [int(p.trace()) for p in (p_in,p_out,p_gamma)]
            row += [nonzero_block(a,lepton[j],anti_sterile,sign*sp.conjugate(ynu)),
                    nonzero_block(delta,quark[i][j],down[i],yd)]
            if np.count_nonzero(a)!=np.count_nonzero(a[np.ix_(lepton[j],anti_sterile)]) or (
                    np.count_nonzero(delta)!=np.count_nonzero(delta[np.ix_(quark[i][j],down[i])])):
                raise ArithmeticError("unexpected additional endpoint matrix elements")
            if row != [3,3,3,3,3]:
                raise ArithmeticError("native charge endpoint ranks changed")
            colour_ranks.append(row)
            # X has the weight of delta = e_i + e_j. On this early
            # weight X vanishes; there is no fibre of the shifted charge.
            # F^2 vanishes on the early Higgs range. Only -2 F X F
            # survives in the projected root-string identity.
            f = vertex.conj().T
            absent = np.array(early_weight)+np.array(late_weight)
            if projector(absent).trace() or np.count_nonzero(f@f@p_in):
                raise ArithmeticError("the two outer root-string terms are not known to vanish")
            if np.count_nonzero(f@p_in-a.conj().T) or np.count_nonzero(f@p_out-delta.conj().T):
                raise ArithmeticError("the lowered endpoint escapes the audited family blocks")
            if np.count_nonzero(f@p_sterile):
                raise ArithmeticError("the audited anti-sterile sector is not annihilated by the lowering field")
            string_checks.append([i+1,j+1])
            if i==0:
                polar_checks.append(_polar_car_check(model,vertex,t3,yukawas))
        ranks.append(colour_ranks)
    if p0.trace()!=6 or n!=96:
        raise ArithmeticError("native probe only applies to the audited 96 fibre")
    return {
        "status":"NATIVE_HIGGS_ENDS_FOUND_CHARGED_MIDDLE_OPEN",
        "native_inputs_sha256":assembler.domain["_sha256"],
        "higgs_builder_sha256":digest(HIGGS),
        "native_gauge_colour_convention":"physical particle fundamental; implemented via the conjugate input in the real lift",
        "rank_per_colour_and_weak_weight":ranks,
        "rank_order":["early_doublet_component","late_quark_component","right_down_component",
                      "neutral_to_early_Higgs","right_down_to_late_Higgs"],
        "neutral_rank":int(p0.trace()),
        "native_higgs_blocks":["anti_nu_R to conjugate L_L: ± conjugate(Ynu)",
                                "d_R to Q_L: Yd"],
        "native_root_string": {
            "affine_identity":"[E_{-alpha,0},[E_{-alpha,0},E_{delta,-1}]]=0",
            "projected_native_residual":"-2 D^* X A^*",
            "family_block_allowed_parameters":0,
            "tested_colour_and_weak_weights":string_checks,
            "reason":"A and D are both full rank; the terms F^2 X and X F^2 vanish on the audited early weight",
            "excluded":"nonzero anti-sterile to Q_L block of delta in the unchanged one-particle fibre with the Higgs block as alpha",
            "gamma_consequence":"both corresponding gamma blocks vanish by [F,X], even if other charged matrices are supplied",
            "extension_test":"if the native fibre reduces the alpha operations and the gauge Cartans, enlarging only homogeneous coloured fields leaves the projected obstruction intact",
        },
        "polar_car_parity":polar_checks,
        "exact_centre_zero_on_native_generators":True,
        "candidate_split":"the colourless weak leg matches Higgs charge blocks; the full colour-changing leg is absent",
        "missing_middle":"charged maps neutral to right-down and early lepton doublet to late quark doublet, with their adjoints and a common native source binding",
        "physical_source_identity_proved":False,
        "original_values":dict.fromkeys(MARKS),
    }


def _polar_car_check(model, vertex, weak_weight, yukawas):
    """Exhaust the parity obstruction for the polar form on this full fibre."""
    import numpy as np
    import sympy as sp
    domain = np.flatnonzero(np.any(vertex,axis=0))
    image = np.flatnonzero(np.any(vertex,axis=1))
    if (len(domain),len(image)) != (48,48) or np.intersect1d(domain,image).size or np.count_nonzero(vertex@vertex):
        raise ArithmeticError("polar Cartan needs a full square-zero source and target split")
    if any(m.rank()!=3 for m in yukawas.values()):
        raise ArithmeticError("polar certificate needs the four full-rank native family matrices")
    families = [[offset+model.IDX[(*species,g)] for g in range(model.NG)]
                for offset in (0,model.NP) for species in sorted({key[:-1] for key in model.IDX})]
    forms = [(name,np.asarray(variant,dtype=complex)) for name,y in yukawas.items()
             for form in (y,y.T,sp.conjugate(y),y.H) for variant in (form,-form)]
    found, used_rows, used_columns = Counter(), set(), set()
    for rows in families:
        for columns in families:
            family = vertex[np.ix_(rows,columns)]
            if not np.any(family):
                continue
            matches = {name for name,original in forms if np.array_equal(family,original)}
            if len(matches)!=1 or used_rows.intersection(rows) or used_columns.intersection(columns):
                raise ArithmeticError("native polar support needs disjoint, verified family blocks")
            found.update(matches)
            used_rows.update(rows)
            used_columns.update(columns)
    if (used_rows!=set(image) or used_columns!=set(domain) or
            found!=Counter({"Yu":6,"Yd":6,"Ye":2,"Ynu":2})):
        raise ArithmeticError("family blocks do not exhaust the whole polar support")
    block = vertex[np.ix_(image,domain)]
    # The exact family ranks give rank 48. Independently check the formula for
    # the Cartan of the polar factor against the original numerical matrices.
    left,_,right = np.linalg.svd(block)
    partial_isometry = np.zeros_like(vertex)
    partial_isometry[np.ix_(image,domain)] = left@right
    metric_cartan = np.zeros(vertex.shape[0])
    metric_cartan[image],metric_cartan[domain]=1,-1
    if np.max(np.abs(partial_isometry@partial_isometry.conj().T-
                     partial_isometry.conj().T@partial_isometry-np.diag(metric_cartan)))>1e-12:
        raise ArithmeticError("native polar Cartan does not match the occupation grading")
    # Any parity-even native CAR word has an even sum of +/-1 adjoint weights.
    return {"weak_weight":weak_weight,"native_polar_rank":48,
            "cartan_eigenvalue_multiplicities":{"-1":48,"1":48},
            "disjoint_full_rank_family_blocks":dict(found),
            "delta_required_adjoint_weight":1,"gamma_required_adjoint_weight":-1,
            "even_car_allowed_odd_weight_space":0,
            "identity":"exp(i*pi*dGamma(H_alpha)) = (-1)^N on the native exterior space",
            "scope":"all parity-even native CAR operators with the canonically polar-normalized alpha; extensions covered only when K stays trivial",
            "constructive_requirement":"a nontrivial -1 sector of K=Pi exp(i*pi H_alpha); an even delta field must exchange the two K sectors"}


def _root_weight(root):
    q = sum(k*x for k,x in zip((-2,-2,-2,3,3,0,0,0),root))//2
    return (q,(root[3]-root[4])//2,(root[0]-root[1])//2,(root[1]-root[2])//2)


def mirror_spectator(root):
    """Test a time-neutral odd factor against the original separated DET16 block."""
    import sympy as sp

    receipt = r104(root)
    if receipt["status"]!="EXACT_R104_FACTOR":
        raise ArithmeticError("the original DET16 block is required for the mirror parity test")
    modes = receipt["modes"]
    phase = sp.exp(sp.I*sp.pi/7)
    vacuum = {0:1/sp.sqrt(2),(1<<modes)-1:phase/sp.sqrt(2)}
    attenuation = sp.Symbol("a", positive=True)
    fixed = thermal(vacuum,vacuum,attenuation)
    if any(sp.simplify(vacuum[state]-fixed[state])!=0 for state in vacuum):
        raise ArithmeticError("DET16 heat fails to fix the coherent ground state")
    for mode in range(modes):
        majorana = {**car_word(vacuum,[(mode,True)]),**car_word(vacuum,[(mode,False)])}
        image = thermal(majorana,vacuum,attenuation)
        if (not all(state.bit_count()%2 for state in majorana) or
                sp.simplify(inner(majorana,image)-attenuation)!=0):
            raise ArithmeticError("original DET16 odd Majorana response changed")

    # All odd states have energy 2 under the pinned R104 Hamiltonian.
    # The only even state with a different energy is its one-dimensional ground space.
    odd_rank = 1 << (modes-1)
    thermal_trace = det_word(modes,"1/7",[[1,"create"],[1,"annihilate"]],[0,"1/3","2/3"])
    neutral_trace = det_word(modes,"1/7",[],[1])
    expected_word = (odd_rank-1)*sp.exp(-2)+(sp.exp(-sp.Rational(4,3))+sp.exp(-sp.Rational(2,3)))/2
    expected_empty = (2*odd_rank-1)*sp.exp(-2)+1
    if (sp.simplify(sp.sympify(thermal_trace["value"])-expected_word)!=0 or
            sp.simplify(sp.sympify(neutral_trace["value"])-expected_empty)!=0):
        raise ArithmeticError("ordered DET16 thermal word lost its independent time weights")
    return {"status":"NO_TIME_NEUTRAL_UNITARY_ODD_SPECTATOR_ON_DET16",
            "ground_energy":"0","odd_energy":receipt["cap"]["temporal_gap"],
            "odd_ground_rank":0,"odd_excited_rank":odd_rank,"even_excited_rank":odd_rank-1,
            "time_commuting_odd_space_dimension":2*odd_rank*(odd_rank-1),
            "maximum_time_commuting_odd_operator_rank":2*(odd_rank-1),
            "odd_ground_return":"exp(-2*t) * ||B phi||^2",
            "unit_norm_odd_leg_required_rest_energy_shift":"-2 for a grade-preserving transition out of the product ground in the DET16 factor",
            "unprojected_thermal_word_control":thermal_trace,
            "unprojected_thermal_empty_control":neutral_trace,
            "checked_single_mode_majorana_ground_returns":modes,
            "scope":"original separated DET16 block with alpha and the time-neutral insertion acting as a spectator; fixed total grade can change the conclusion"}


def replay(root, domain):
    root, domain = root.resolve(), domain.resolve()
    # Isolate the old source's generic module names from the search engine.
    child = """
from pathlib import Path
import json,sys
for name in sys.argv[1:5]:
    sys.path.insert(0,name)
from native_weak_bridge import analyze
print(json.dumps(analyze(Path(sys.argv[5]))))
"""
    paths = [root/RUNTIME,root/COMPARATOR,root/"verification",Path(__file__).parent]
    run = subprocess.run([sys.executable,"-B","-c",child,*map(str,paths),str(domain)],
                         cwd=root,capture_output=True,text=True,timeout=60)
    if run.returncode:
        raise RuntimeError(run.stderr.strip().splitlines()[-1] if run.stderr else "native weak check failed")
    return {**json.loads(run.stdout),"det16_spectator":mirror_spectator(root),
            "source_refs":{str(p):digest(root/p) for p in INPUTS}}
