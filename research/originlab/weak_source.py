"""Gauge-resolved affine benchmark and same-source R113 word diagnostics."""

from collections import Counter
from dataclasses import dataclass
from itertools import product

import sympy as sp

from affine_operators import (LatticeFock, MARKED_GAMMA, MARKED_P1, ROOTS, ROOT_SET, WEAK_ROOT,
                              add, charge_signature, closure, dot, integral, lattice_frame, neg)
from reconstruction import matrix, require, scalar


def axis(index):
    return tuple(2 if i == index else 0 for i in range(8))


@dataclass(frozen=True)
class GaugeChannel:
    colour: str
    weak: int
    charge: int  # 6Y

    def casimirs(self):
        return (sp.Rational(4, 3) if self.colour != "1" else sp.S.Zero,
                sp.Rational(3, 4) if self.weak == 2 else sp.S.Zero)

    def centre_condition(self):
        triality = {"1": 0, "3": 1, "bar3": 2}[self.colour]
        return (self.charge + 3*(self.weak-1) + 2*triality) % 6 == 0


# Intervals are in chronological order, from the right end toward the left.
WEAK_CHANNELS = (GaugeChannel("1", 2, 3), GaugeChannel("3", 2, 1))
CHANNELS = {
    "spinor": (GaugeChannel("bar3", 1, 2), GaugeChannel("1", 1, 0)),
    "d8": WEAK_CHANNELS,
    "weak_spinor_minus": WEAK_CHANNELS,
    "weak_spinor_plus": WEAK_CHANNELS,
}
NETWORKS = ("identity", "weak_first", "weak_inverse", "colour_first", "colour_second",
            "noncommuting_weak_pair", "colour_hypercharge")
PROBES = ("early_colour", "early_weak", "late_colour", "late_weak")
WINDING_PROBES = ("early_half", "late_half", "both_halves")
COLOUR_SIMPLE_ROOTS = (add(axis(0), neg(axis(1))), add(axis(1), neg(axis(2))))
COWEIGHT = neg(add(axis(7), axis(7)))
WEAK_TORUS_ENTRY = (3+4*sp.I)/5


def triangle_components():
    """Search all roots for a weak doublet leg with the middle charge held fixed."""
    tensors = {name: [] for name in CHANNELS}
    colour_steps = (*COLOUR_SIMPLE_ROOTS, *(neg(r) for r in COLOUR_SIMPLE_ROOTS))
    candidates = []
    for root in ROOTS:
        _, triality, charge = charge_signature(root)
        if (triality == 0 and 2*charge == 1 and abs(dot(root, WEAK_ROOT)) == 1 and
                all(dot(root, step) == 0 and add(root, step) not in ROOT_SET for step in colour_steps)):
            candidates.append(root)
    for i in range(3):
        middle = add(axis(i), axis(5))
        tensors["spinor"].append((neg(MARKED_P1), middle, add(MARKED_P1, neg(middle))))
        for right in candidates:
            total = add(middle, right)
            if total not in ROOT_SET:
                continue
            if all(entry % 2 == 0 for entry in right):
                name = "d8"
            else:
                winding = dot(COWEIGHT, right)
                require(abs(winding) == 1, "unexpected spinor winding in the fixed-middle search")
                name = "weak_spinor_minus" if winding < 0 else "weak_spinor_plus"
            tensors[name].append((neg(total), middle, right))
    require([len(terms) for terms in tensors.values()] == [3, 6, 6, 6],
            "fixed-middle colour-singlet weak legs did not fill the expected gauge orbits")
    return tensors


def native_d8_legs():
    # The Higgs quaternion chart lists the + weak weight first. Keep its
    # ordering explicit and independent of the order of root enumeration.
    return sorted(triangle_components()["d8"],
                  key=lambda word:(-word[1][0],-word[1][1],-dot(word[2],WEAK_ROOT)))


def charge_classes(roots):
    counts = Counter()
    for root in roots:
        _, triality, hypercharge = charge_signature(root)
        q = 6*hypercharge
        require(q.denominator == 1 and q % 2 == dot(root, WEAK_ROOT) % 2 and q % 3 == triality,
                "incompatible centre character in the fixed D8 dictionary")
        counts[int(q) % 6] += 1
    return dict(sorted(counts.items()))


def winding_characters(words):
    """Use the original relative coweight -2e8, not the SM U(1) charge."""
    labels = {(dot(COWEIGHT, word[2]), dot(COWEIGHT, neg(word[0]))) for word in words}
    require(len(labels)==1, "coweight must be fixed on each complete gauge orbit")
    early,late = labels.pop()
    return (early,late),(-1)**(early%2),(-1)**(late%2),(-1)**((early+late)%2)


def _joint_early_phase(root, power):
    return (WEAK_TORUS_ENTRY**(power*dot(root, WEAK_ROOT))*
            sp.I**((power*dot(COWEIGHT, root)) % 4))


def joint_early_character(words, power=1):
    """One diagonal weak and Cartan insertion, in the first interval only."""
    # Both actions commute with the full Cartan. The normalized coefficients
    # of every component are equal, as checked with the mode operators below.
    return sp.simplify(sum(_joint_early_phase(word[2], power) for word in words)/len(words))


def joint_response_matrix(tensors):
    return sp.Matrix([[joint_early_character(words,power) for words in tensors.values()]
                      for power in range(5)])


def _radial_source_sample(x, y, held_out=False):
    # The x interval carries degree m; the y interval carries degree n.
    # Sum min(m,n) x**m y**n with the cocycle corrected and no clock phases.
    kernel = x*y/((1-x)*(1-y)*(1-x*y))
    return {"left_radial_factor":str(x),"right_radial_factor":str(y),
            "affine_bare_response":str(kernel),"held_out":held_out}


def _commutator_action(tensor, generator, fock, frame):
    result = Counter()
    v = integral(frame, generator)
    for word, coefficient in tensor.items():
        for i, root in enumerate(word):
            target = add(generator, root)
            # The diagonal Cartan branch never occurs for the six gauge step roots.
            require(root != neg(generator), "unexpected Cartan term in triangle action")
            if target in ROOT_SET:
                image = (*word[:i], target, *word[i+1:])
                result[image] += coefficient*fock.sign(v, integral(frame, root))
    return {word: value for word, value in result.items() if value}


def _unitary_link(data):
    w, c, u = matrix(data["weak"]), matrix(data["colour"]), scalar(data["u1"])
    for group, size in ((w, 2), (c, 3)):
        require(group.shape == (size, size) and group.H*group == sp.eye(size) and group.det() == 1,
                "weak and colour links must lie in SU(2) and SU(3)")
    require(sp.expand(u*sp.conjugate(u)) == 1, "hypercharge link must have modulus one")
    return w, c, u


def transport_character(name, first, second):
    """Response to gauge actions on the two intermediate affine channels."""
    require(name in CHANNELS, "unknown gauge channel path")
    w1,c1,u1 = first
    w2,c2,u2 = second
    if name == "spinor":
        return sp.simplify(sp.conjugate(sp.trace(c1))*u1**2/3)
    return sp.simplify(sp.trace(c2)*sp.trace(w2*w1)*u1**3*u2/6)


def example_networks():
    weak, colour = sp.eye(2), sp.eye(3)
    identity = (weak, colour, sp.S.One)
    rotation = sp.diag(sp.I, -sp.I)
    transverse = sp.Matrix([[0, sp.I], [sp.I, 0]])
    colour_flip = sp.diag(-1, -1, 1)
    return {
        "identity": (identity, identity),
        "weak_first": ((rotation, colour, 1), identity),
        "weak_inverse": ((rotation, colour, 1), (rotation.H, colour, 1)),
        "colour_first": ((weak, colour_flip, 1), identity),
        "colour_second": (identity, (weak, colour_flip, 1)),
        "noncommuting_weak_pair": ((rotation, colour, 1), (transverse, colour, 1)),
        "colour_hypercharge": ((weak, colour_flip, sp.I), (weak, colour_flip, -sp.I)),
    }


def _group_casimir(fock, frame, state):
    # The lattice cocycle gives E_r(0)^* = -E_{-r}(0). For A2 the
    # Cartan term is half the norm of the colour weight projection.
    require(bool(state) and len({key[0] for key in state}) == 1, "Casimir needs one nonzero weight")
    charge = next(iter(state))[0]
    h = fock.pair(integral(frame, WEAK_ROOT), charge)
    colour_roots = [add(axis(i), neg(axis(j))) for i,j in ((0,1), (1,2), (0,2))]
    differences = [fock.pair(integral(frame, root), charge) for root in colour_roots]

    def action(roots, cartan):
        steps = (integral(frame,r) for root in roots for r in (root,neg(root)))
        return fock.combine((cartan,state), *((sp.Rational(-1,2),
                            fock.vertex(v,0,fock.vertex(neg(v),0,state))) for v in steps))

    return (action(colour_roots, sum(sp.Rational(k*k,6) for k in differences)),
            action([WEAK_ROOT], sp.Rational(h*h,4)))


def _native_d8_string(fock, frame, components):
    checks = 0
    for left,middle,right in components:
        lower, delta, gamma = [integral(frame,root) for root in (neg(right),neg(left),middle)]
        sign = fock.sign(lower,delta)
        vacuum = fock.vacuum()
        for state in (vacuum,fock.vertex(neg(lower),-1,vacuum)):
            def f(vector):
                return fock.vertex(lower,0,vector)
            def x(vector):
                return fock.vertex(delta,-1,vector)
            y = fock.combine((1,f(x(state))),(-1,x(f(state))))
            require(y == fock.combine((sign,fock.vertex(gamma,-1,state))),
                    "first D8 root-string commutator changed")
            require(not fock.combine((1,f(f(x(state)))),(-2,f(x(f(state)))),(1,x(f(f(state))))),
                    "second D8 root-string bracket must vanish")
            checks += 2
    return checks


def parity_completion():
    """The minimal horizontal escape when the alpha doublet is retained."""
    identity = sp.eye(3)
    alpha = identity[:,2]*identity[1,:]
    gamma = identity[:,1]*identity[0,:]
    delta = identity[:,2]*identity[0,:]
    h = alpha*alpha.H-alpha.H*alpha
    p = sp.diag(-1,-1,-1)
    centre = sp.diag(1,-1,-1)
    mismatch = p*centre
    require(h == sp.diag(0,-1,1) and alpha*gamma-gamma*alpha == delta and
            all(e*e.H-e.H*e==diag and diag*e-e*diag==2*e
                for e,diag in ((alpha,h),(gamma,sp.diag(-1,1,0)),(delta,sp.diag(-1,0,1)))) and
            alpha.H*delta-delta*alpha.H == gamma and delta.H*alpha-alpha*delta.H == gamma.H and
            mismatch*alpha == alpha*mismatch and mismatch*delta == -delta*mismatch and
            delta.H*gamma*alpha == sp.zeros(3) and delta.H*alpha*gamma == identity[:,0]*identity[0,:],
            "horizontal A2 extension failed")
    # The affine NS vector realization uses a pair for delta. Reuse the
    # complete CAR construction already used for the little coherence tests.
    from coherence import car_word

    def operator(word):
        out = sp.zeros(8)
        for occupation in range(8):
            for target,value in car_word({occupation:1},word).items():
                out[target,occupation] = value
        return out

    alpha_car = operator([(1,True),(2,False)])
    gamma_car = operator([(0,True),(2,True)])
    delta_car = alpha_car*gamma_car-gamma_car*alpha_car
    cartan = alpha_car*alpha_car.H-alpha_car.H*alpha_car
    omega = sp.eye(8)[:,0]
    early, late = gamma_car*omega, delta_car*omega
    parity = sp.diag(*((-1)**i.bit_count() for i in range(8)))
    report = checked_parity_sector({
        "parity_form":parity,"alpha_cartan_form":cartan,"alpha_raising_form":alpha_car,
        "alpha_lowering_form":alpha_car.H,"delta_form":delta_car,"delta_adjoint_form":delta_car.H,
        "neutral_projector_form":omega*omega.H,"early_projector_form":early*early.H,
        "late_projector_form":late*late.H,
    },sp.eye(8),lambda value: value)
    require(report["status"]=="NONZERO_SUPPLIED_PARITY_BRIDGE" and
            report["alpha_then_gamma_return_rank"]==0 and report["gamma_then_alpha_return_rank"]==1 and
            (omega.H*delta_car.H*alpha_car*gamma_car*omega)[0]==1,
            "vacuum-supported D8 pair realization changed")
    return {"status":"EXACT_ONE_EXTRA_ALPHA_SINGLET_CONTROL",
            "source_doublet_dimension":2,"extra_alpha_singlet_dimension":1,
            "lie_actions":"alpha: lower to upper doublet weight; gamma: singlet to lower; delta: singlet to upper",
            "K_eigenvalues":[-1,1,1],"intermediate_probe_order_coefficients":["0","1"],
            "vacuum_pair_order_coefficients":["0","1"],"vacuum_pair_endpoint_test":report,
            "scope":"one horizontal A2 block and a conditional spectator; source, affine central term and locality are open"}


def checked_parity_sector(data, metric, read):
    """Check for the missing even, odd-Cartan field before outer selection."""
    identity = sp.eye(metric.rows)
    parity, h, raising, lowering, delta, delta_star, neutral, early, late = [read(data[key]) for key in
        ("parity_form","alpha_cartan_form","alpha_raising_form","alpha_lowering_form","delta_form",
         "delta_adjoint_form","neutral_projector_form","early_projector_form","late_projector_form")]
    metric_inv = metric.inv()

    def adjoint(op):
        return metric_inv*op.H*metric

    require(parity==adjoint(parity) and parity**2==identity and
            h==adjoint(h) and parity*h==h*parity, "joint parity and alpha Cartan action is invalid")
    require(lowering==adjoint(raising) and delta_star==adjoint(delta), "original field adjoints do not match")
    levels = h.eigenvals()
    require(all(value.is_Integer for value in levels), "compact alpha Cartan requires integer weights")
    # Positive metric gives a diagonalizable h. Polynomial evaluation
    # implements exp(i*pi*h) exactly, including degeneracies.
    centre = sp.zeros(metric.rows)
    for value in levels:
        projection = identity
        for other in levels:
            if other!=value:
                projection = projection*(h-other*identity)/(value-other)
        centre += (-1)**(value%2)*projection
    require(all(p==adjoint(p) and p*p==p and p*parity==parity*p and p*h==h*p
                for p in (neutral,early,late)), "endpoint projections must reduce parity and the alpha Cartan")
    k = parity*centre
    gamma = lowering*delta-delta*lowering
    residuals = {
        "alpha_parity":(parity*raising-raising*parity).rank(),
        "alpha_cartan":(h*raising-raising*h-2*raising).rank(),
        "alpha_commutator":(raising*lowering-lowering*raising-h).rank(),
        "delta_parity":(parity*delta-delta*parity).rank(),
        "delta_cartan":(h*delta-delta*h-delta).rank(),
        "root_string_raise":(raising*delta-delta*raising).rank(),
        "root_string_lower":(lowering*gamma-gamma*lowering).rank(),
        "root_string_return":(raising*gamma-gamma*raising-delta).rank(),
        "mismatch_sector":(k*delta+delta*k).rank(),
    }
    block = late*delta*neutral
    preparations = early*raising*neutral
    early_path = late*gamma*preparations
    reversed_path = late*raising*gamma*neutral
    status = "NO_EVEN_ODD_WEIGHT_ENDPOINT"
    if any(residuals.values()):
        status = "INCOMPATIBLE_SUPPLIED_ALPHA_WORDS"
    elif block.rank():
        status = "NONZERO_SUPPLIED_PARITY_BRIDGE"
    return {"status":status,"k_minus_rank":(identity-k).rank(),"delta_endpoint_rank":block.rank(),
            "early_k_minus_rank":(early*(identity-k)).rank(),"residual_ranks":residuals,
            "alpha_then_gamma_return_rank":(adjoint(block)*early_path).rank(),
            "gamma_then_alpha_return_rank":(adjoint(block)*reversed_path).rank(),
            "physical_source_identity_proved":False}


def affine_controls(mode_cutoff=6):
    """Replay the marked extension and every fixed-middle weak orbit."""
    require(type(mode_cutoff) is int and 1 <= mode_cutoff <= 6, "mode cutoff must be 1 through 6")
    frame,_ = lattice_frame()
    fock = LatticeFock([[dot(a,b) for b in frame] for a in frame])
    vacuum = fock.vacuum()
    tensors = triangle_components()
    even = {r for r in ROOTS if dot(r,WEAK_ROOT) % 2 == 0}
    odd = ROOT_SET-even
    alpha = tensors["d8"][0][2]
    partners = {frozenset((r, add(r, neg(WEAK_ROOT)))) for r in odd if dot(r,WEAK_ROOT) == 1}
    require(len(odd) == 112 and len(partners) == 56 and set().union(*partners) == odd and
            closure((*even, alpha)) == ROOT_SET and all(len(pair) == 2 for pair in partners),
            "weak doublet or closure census changed")
    require(all((sp.Rational(dot(r,WEAK_ROOT)**2,4), sp.Rational(1)-sp.Rational(dot(r,WEAK_ROOT)**2,4))
                == (sp.Rational(1,4),sp.Rational(3,4)) for r in odd), "wrong affine weak factor weights")
    gauges = (*COLOUR_SIMPLE_ROOTS, WEAK_ROOT)
    invariance_count = mode_count = casimir_count = joint_count = 0
    component_phases = {}
    winding = {}
    signatures = {}
    readout = joint_response_matrix(tensors)
    for column, (name, components) in enumerate(tensors.items()):
        terms = {}
        joint = readout[1,column]
        fixed_middle = [word for word in components if word[1] == MARKED_GAMMA]
        glue = {(charge_signature(word[2])[0],charge_signature(neg(word[0]))[0]) for word in components}
        require(len(glue) == 1, "glue grading must be constant on the gauge orbit")
        signatures[name] = {"fixed_middle_right_roots_doubled": [list(word[2]) for word in fixed_middle],
                            "interval_glue": list(next(iter(glue))), "joint_early_ratio": str(joint)}
        joint_amplitudes = {pair:sp.zeros(readout.rows,1)
                            for pair in product(range(1,mode_cutoff+1),repeat=2)}
        for word in components:
            require(add(add(*word[:2]), word[2]) == (0,)*8 and
                    all(dot(a,b) == -1 for a,b in ((word[0],word[1]),(word[0],word[2]),(word[1],word[2]))),
                    "neutral A2 triangle required")
            left, middle, right = [integral(frame,r) for r in word]
            phase = fock.sign(middle,right)*fock.sign(left,neg(left))
            terms[word] = phase
            insertions = sp.Matrix([sp.simplify(_joint_early_phase(word[2],power))
                                    for power in range(readout.rows)])
            for m,n in joint_amplitudes:
                early = fock.vertex(right,-n,vacuum)
                late = fock.vertex(middle,n-m,early)
                out = fock.vertex(left,m,late)
                bare = phase*fock.inner(vacuum,out)
                require(bare == min(m,n), "cocycle-normalized mode disagreement")
                mode_count += 1
                joint_amplitudes[m,n] += insertions*bare
                # Resolve intermediate irreps by inserting the Casimir.
                # After the middle leg the left cap selects a channel of
                # the full product, which may contain other representations.
                for early_c, late_c, early_expected, late_expected in zip(
                    _group_casimir(fock,frame,early), _group_casimir(fock,frame,late),
                    CHANNELS[name][0].casimirs(), CHANNELS[name][1].casimirs()
                ):
                    early_word = fock.vertex(left,m,fock.vertex(middle,n-m,early_c))
                    late_word = fock.vertex(left,m,late_c)
                    require(phase*fock.inner(vacuum,early_word) == early_expected*min(m,n) and
                            phase*fock.inner(vacuum,late_word) == late_expected*min(m,n),
                            "interval Casimir disagrees with the actual affine mode")
                    casimir_count += 2
        for (m,n), amplitude in joint_amplitudes.items():
            require((amplitude/len(components)-readout[:,column]*min(m,n)).applyfunc(sp.simplify)
                    == sp.zeros(readout.rows,1),
                    "joint weak and coweight action changed the affine mode response")
            joint_count += readout.rows
        for gauge in gauges:
            for direction in (gauge,neg(gauge)):
                require(not _commutator_action(terms,direction,fock,frame), "not a global gauge singlet")
                invariance_count += 1
        component_phases[name] = list(terms.values())
        grades,*phase = winding_characters(components)
        winding[name] = {"interval_winding":list(grades),"inserted_half_turns":dict(zip(WINDING_PROBES,phase))}
        require(all(ch.centre_condition() for ch in CHANNELS[name]), "interval transport fails the Z6 quotient")
    require([winding[name]["interval_winding"] for name in CHANNELS] == [[-1,-1],[0,0],[-1,-1],[1,1]]
            and [entry["interval_glue"] for entry in signatures.values()] == [[1,3],[2,0],[3,1],[3,1]],
            "fixed-middle winding or glue action changed")
    def exchange_last_slots(root):
        return (*root[:6], root[7], root[6])

    require({tuple(map(exchange_last_slots,word)) for word in tensors["weak_spinor_minus"]} ==
            set(tensors["weak_spinor_plus"]), "weak spinor orbits are not related by the e7-e8 reflection")
    determinant = sp.simplify(readout[:4,:].det())
    require(determinant != 0 and len(set(readout[1,:])) == 4,
            "the source probe cannot distinguish the four fixed-middle paths")
    widths = (sp.Rational(1,2),sp.Rational(1,3))
    training_grid = [_radial_source_sample(x,y) for x,y in product(widths,repeat=2)]
    holdout = _radial_source_sample(sp.Rational(1,4),sp.Rational(1,5),held_out=True)

    networks = example_networks()
    characters = {k:{name:str(transport_character(name,*pair)) for name in CHANNELS}
                  for k,pair in networks.items()}
    require([list(characters[key].values()) for key in NETWORKS[:5]] ==
            [["1"]*4,["1","0","0","0"],["1"]*4,["-1/3","1","1","1"],["1","-1/3","-1/3","-1/3"]],
            "independent interval character controls did not separate the words")
    require(list(characters["noncommuting_weak_pair"].values()) == ["1","0","0","0"] and
            set(characters["colour_hypercharge"].values()) == {"1/3"}, "combined link controls changed")
    return {
        "status": "EXACT_AFFINE_GAUGE_CHANNELS",
        "root_centre_counts": charge_classes(ROOTS), "weak_doublets": len(partners),
        "affine_weights": ["1/4","3/4","1"], "affine_twists": ["i","-i","1"],
        "colour_and_weak_invariance_checks": invariance_count, "mode_checks": mode_count,
        "casimir_checks": casimir_count, "component_phases": component_phases,
        "d8_root_string_checks": _native_d8_string(fock,frame,tensors["d8"]),
        "intervals": {name:[{"colour":q.colour,"weak_dimension":q.weak,"six_Y":q.charge,
                             "casimirs":list(map(str,q.casimirs()))} for q in pair]
                      for name,pair in CHANNELS.items()},
        "character_controls": characters,
        "original_coweight_controls":winding,
        "fixed_middle_search":{"examined_roots":len(ROOTS),"paths":signatures,
                               "weak_and_coweight_mode_checks":joint_count,
                               "column_order":list(tensors),"joint_power_response":
                                   [[str(value) for value in row] for row in readout.tolist()],
                               "powers_zero_to_three_determinant":str(determinant)},
        "polar_parity_necessary_extension_control":parity_completion(),
        "next_original_observation": {
            "status":"REQUIRES_ORIGINAL_SOURCE_BINDING",
            "experiment":"insert the weak torus element and quarter coweight on the early interval, with the same fields and native dressing",
            "first_path_to_compare_with_the_marked_word":"weak_spinor_minus",
            "weak_torus_diagonal":[str(WEAK_TORUS_ENTRY),str(sp.conjugate(WEAK_TORUS_ENTRY))],
            "quarter_coweight":"exp(i*pi*H_(-2e8)/2)","chronological_intervals":["y","x"],
            "word_forms":"absolute preselection forms of insertion powers zero through four in a common word basis",
            "sample_grid":[*training_grid,holdout],
            "source_requirements":["native ordered weighted form","physical to radial time and fixed gauge frames",
                                   "one outer selector","original absolute empty reference"],
        },
        "physical_source_identity_proved": False,
    }


def checked_responses(data, contract):
    """Compare gauge responses using inserted actions on one supplied source.

    The affine factors are predictions for external gauge actions in the
    bound frames. A dynamical link correlated with the rest of the field
    history need not factor out of the original measure.
    """
    samples = data["samples"]
    require(isinstance(samples,list) and samples, "use a grid and an independently reserved sample")
    seen = set()
    training, held_out = set(), set()
    report = []
    networks = example_networks()
    triangles = triangle_components()
    paths = data.get("comparison_paths", ["spinor","d8"])
    require(isinstance(paths,list) and paths and paths[0]=="spinor" and 2<=len(paths)<=len(CHANNELS)
            and all(isinstance(name,str) and name in CHANNELS for name in paths)
            and len(set(paths))==len(paths), "select the marked word and at least one distinct fixed-middle path")
    has_joint = any("joint_early_forms" in sample for sample in samples)
    if has_joint:
        require(isinstance(data.get("joint_early_insertion_binding"),str) and
                data["joint_early_insertion_binding"].strip() and
                all("joint_early_forms" in sample for sample in samples),
                "bind the same original early Cartan and weak insertion on every sample")
        response = joint_response_matrix(triangles)
        inverse = response[:4,:].inv()
    for sample in samples:
        x,y = (scalar(sample[axis]) for axis in ("left_radial_factor","right_radial_factor"))
        holdout = sample["held_out"]
        require(type(holdout) is bool and x.is_Rational and y.is_Rational and 0<x<1 and 0<y<1
                and (x,y) not in seen, "use unique ordered positive intervals")
        seen.add((x,y))
        (held_out if holdout else training).add((x,y))
        row = {"radial_factors":[str(x),str(y)],"held_out":holdout}
        if has_joint:
            forms = sample["joint_early_forms"]
            require(isinstance(forms,list) and len(forms)==5,
                    "supply absolute forms for insertion powers zero through four")
            measured = sp.Matrix([contract(form) for form in forms])
            weights = (inverse*measured[:4,:]).applyfunc(sp.simplify)
            row["joint_early"] = {"absolute_power_words":list(map(str,measured)),
                                 "absolute_path_coefficients":dict(zip(triangles,map(str,weights))),
                                 "power_four_residual":str(sp.simplify(measured[4]-(response[4,:]*weights)[0]))}
        require(set(sample["words"])==set(paths), "each sample must supply the chosen fixed-middle paths")
        for name in paths:
            observed = sample["words"][name]
            channels = CHANNELS[name]
            bare = contract(observed["bare_form"])
            targets = [value for channel in channels for value in channel.casimirs()]
            casimirs = {p:contract(observed["casimir_forms"][p]) for p in PROBES}
            require(set(observed["twisted_words"]) == set(NETWORKS)-{"identity"},
                    "supply early, late, cancelling and held-out action probes")
            twists = {}
            for key, values in observed["twisted_words"].items():
                supplied = tuple(_unitary_link(link) for link in values["links"])
                require(supplied == networks[key], "link actions must match the two specified group elements in order")
                prediction = transport_character(name,*supplied)
                measured = contract(values["word_form"])
                twists[key] = {"absolute":str(measured),"expected_ratio":str(prediction),
                               "residual":str(sp.simplify(measured-bare*prediction))}
            winding = {}
            if "winding_forms" in observed:
                require(isinstance(data.get("original_coweight_binding"),str) and data["original_coweight_binding"].strip()
                        and set(observed["winding_forms"])==set(WINDING_PROBES),
                        "bind the same original coweight and supply its early, late and paired forms")
                _,*phases = winding_characters(triangles[name])
                for key,phase in zip(WINDING_PROBES,phases):
                    value = contract(observed["winding_forms"][key])
                    winding[key] = {"absolute":str(value),"expected_ratio":str(phase),
                                    "residual":str(sp.simplify(value-bare*phase))}
            residuals = {p:str(sp.simplify(casimirs[p]-bare*target)) for p,target in zip(PROBES, targets)}
            row[name] = {"absolute_bare":str(bare), "degenerate":bare==0,
                         "absolute_casimir_words":{p:str(v) for p,v in casimirs.items()},
                         "casimir_residuals":residuals,"transport":twists,"winding":winding}
        report.append(row)
    x_values,y_values = ({p[i] for p in training} for i in (0,1))
    require(len(x_values) >= 2 and len(y_values) >= 2 and training == set(product(x_values,y_values))
            and any(x not in x_values and y not in y_values for x,y in held_out),
            "supply a complete calibration rectangle and a holdout with both widths unseen")
    checks = [row[name] for row in report for name in paths]
    if any(word["winding"] for word in checks):
        require(all(word["winding"] for word in checks), "winding comparison needs every path on the whole grid")
    mismatch = any(any(residual != "0" for residual in word["casimir_residuals"].values()) or
                   any(probe["residual"] != "0" for values in (word["transport"],word["winding"])
                       for probe in values.values()) for word in checks)
    joint_mismatch = has_joint and any(row["joint_early"]["power_four_residual"] != "0" for row in report)
    status = "MATCHES_SUPPLIED_AFFINE_RESPONSES"
    if mismatch or joint_mismatch:
        status = "DISTINGUISHED_SUPPLIED_GAUGE_RESPONSES"
    elif any(check["degenerate"] for check in checks):
        status = "INCONCLUSIVE_ZERO_BARE_WORD"
    return {"status":status,"samples":report,"comparison_paths":paths,"physical_source_identity_proved":False,
            "scope":"affine response comparison on the supplied whole words with externally inserted group actions"}
