"""Generate and check the coherent two-copy source family from the F1 report."""

import json
from itertools import product
from pathlib import Path
from tempfile import TemporaryDirectory

import sympy as sp

from reconstruction import require
from source import COMPONENTS, evaluate
from verify import digest


def gaussian_covariance():
    """Pure centered Gaussian branch in the report's Fourier convention."""
    a,b,c=sp.symbols("a b c",real=True)
    covariance=sp.Matrix([[a,c],[c,b]])
    rotation=sp.Matrix([[0,1],[-1,0]])
    constraints=list(rotation*covariance*rotation.T-covariance)
    solutions=sp.solve([*constraints,covariance.det()-1/(16*sp.pi**2)],(a,b,c))
    variance,_,_=next(row for row in solutions if row[0].is_positive)
    identity=sp.eye(2)
    mixer=sp.BlockMatrix([[identity,identity],[identity,-identity]]).as_explicit()/sp.sqrt(2)
    joint=sp.eye(4)*variance
    require(mixer*joint*mixer.T == joint, "two-copy Gaussian covariance was not invariant")
    embed=sp.eye(4)-sp.ones(4)/4
    relative=sp.simplify(2*variance*embed/4)
    require(relative==embed/(8*sp.pi), "relative covariance differs from the K4 coefficient")
    return {
        "one_side_variance":str(variance),"raw_difference_variance":str(2*variance),
        "k4_covariance":[[str(v) for v in row] for row in relative.tolist()],
        "mixed_countermodel_variance":str(3*variance),
        "scope":"independent primitive two-copy Gaussian state, purity and fixed Fourier pairing",
    }


def number_sectors(limit):
    """The nonzero beam splitter amplitudes of |n,0> impose energy additivity."""
    failures=[]
    for n in range(limit+1):
        weights=[sp.binomial(n,k)/2**n for k in range(n+1)]
        require(sum(weights)==1 and all(w>0 for w in weights), "incomplete beam splitter sector")
        ratio=sp.Symbol("r",positive=True)
        thermal=lambda j: (1-ratio)*ratio**j
        require(all(sp.expand(thermal(k)*thermal(n-k)-thermal(0)*thermal(n))==0
                    for k in range(n+1)), "geometric weights differ within a fixed number sector")
        for k in range(n+1):
            shifted=lambda m: m+4*m*(m-1)
            if shifted(k)+shifted(n-k)!=shifted(n)+shifted(0):
                failures.append([n,k])
    require(failures[0]==[2,1], "same-quarter nonlinear generator must fail at n=2")
    return {"tested_total_number":limit,"relative_generator":"omega * N",
            "first_nonlinear_alias_failure":failures[0],
            "absolute_vacuum_energy_selected":False,"physical_rate_selected":False}


def decoder_controls():
    """Check exact coherence and additive-energy transfer on the four-code subspace."""
    encoded=sp.zeros(16,4)
    decoder=sp.zeros(16)
    for j,k in product(range(4),repeat=2):
        decoder[4*j+(k-j)%4,4*j+k]=1
        if j==k:
            encoded[4*j+k,j]=1
    decoded=decoder*encoded
    target=sp.kronecker_product(sp.eye(4),sp.eye(4)[:,0])
    require(decoder.H*decoder==sp.eye(16) and decoded==target, "decoder is not a coherent isometry")
    # Conjugating matrix units proves multiplicativity on the full logical algebra.
    x=sp.zeros(4)
    for j in range(4):
        x[(j+1)%4,j]=1
    z=sp.diag(1,sp.I,-1,-sp.I)
    lift=lambda matrix: decoder.H*sp.kronecker_product(matrix,sp.eye(4))*decoder
    require(lift(x*z)==lift(x)*lift(z) and lift(x)*lift(z)!=lift(z)*lift(x),
            "decoder lost a logical product")
    n=sp.diag(*range(4))
    additive=sp.kronecker_product(n,sp.eye(4))+sp.kronecker_product(sp.eye(4),n)
    difference=encoded.H*additive*encoded-decoded.H*additive*decoded
    require(difference==n and difference.rank()==3, "raw decoder energy deficit")
    q=sp.Symbol("q",positive=True)
    total_transfer=sp.diag(*(q**(j+k) for j,k in product(range(4),repeat=2)))
    before=encoded.H*total_transfer*encoded
    raw_after=decoded.H*total_transfer*decoded
    doubled_after=(decoded.H*total_transfer**2*decoded)
    require(before!=raw_after and before==doubled_after, "decoder transfer scale control")
    # Conservation with a vacuum reservoir forces orthogonal records for
    # different logical numbers. Tracing the reservoir removes coherence.
    reservoir=sp.eye(4)
    initial=sp.kronecker_product(encoded,reservoir[:,0])
    conserved=sp.Matrix.hstack(*(sp.kronecker_product(decoded[:,j],reservoir[:,j]) for j in range(4)))
    number_with_reservoir=sp.kronecker_product(additive,reservoir)+sp.kronecker_product(sp.eye(16),n)
    require(conserved.H*conserved==sp.eye(4)
            and initial.H*number_with_reservoir*initial==conserved.H*number_with_reservoir*conserved,
            "vacuum reservoir energy balance")
    for j,k in product(range(4),repeat=2):
        pair_operator=conserved[:,j]*conserved[:,k].H
        reduced=sp.Matrix(4,4,lambda a,b: sum(pair_operator[16*a+m,16*b+m] for m in range(16)))
        require(reduced==(reservoir[:,j]*reservoir[:,k].H if j==k else sp.zeros(4)),
                "tracing reservoir did not produce the predicted logical channel")
    return {"logical_product_preserved":True,"bare_energy_deficit":[0,1,2,3],
            "vacuum_reservoir_energy_deficit":[0,0,0,0],
            "off_diagonal_coherence_after_tracing_reservoir":"0",
            "bare_transfer_intertwines":False,"rescaled_output_transfer_intertwines":True,
            "required_output_step":"q^2",
            "energy_and_coherence_require_common_dressed_observables":True}


def marked_operators():
    swap=sp.zeros(4)
    for j in range(4):
        swap[j^1,j]=1
    parity=sp.diag(1,-1,1,-1)
    phase=sp.I*swap*parity
    require(swap**2==phase**2==sp.eye(4) and swap*phase==sp.I*parity
            and phase*swap==-sp.I*parity, "coherent ordered word algebra")
    return swap,phase,parity


def words_control():
    a,b,parity=marked_operators()
    vacuum=sp.eye(4)[:,0]
    cap=vacuum*vacuum.H
    alphabet=(a,b)
    # The A/B word orbit lies in one even/odd pair. Compare with its exact
    # two-direction restriction, leaving the phase of every word intact.
    checked=0
    for length in range(7):
        for word in product(range(2),repeat=length):
            full,small=sp.eye(4),sp.eye(2)
            for letter in word:
                full*=alphabet[letter]
                small*=alphabet[letter][:2,:2]
            require(full[0,0]==small[0,0],"tail truncation changed an ordered word")
            checked+=1
    require(cap*a*cap==cap*b*cap==sp.zeros(4) and cap*a*b*cap==sp.I*cap,
            "intermediate cap must fail the ordered word test")
    # Two even Majorana bilinears on a three-mode exterior space have the
    # same vacuum response as the two-dimensional A/B representation.
    z=sp.diag(1,-1)
    x=sp.Matrix([[0,1],[1,0]])
    y=sp.Matrix([[0,-sp.I],[sp.I,0]])
    gammas=[sp.kronecker_product(*([z]*j+[axis]+[sp.eye(2)]*(2-j)))
            for j in range(3) for axis in (x,y)]
    even_a=sp.I*gammas[0]*gammas[2]
    even_b=sp.I*gammas[1]*gammas[2]
    fermion_parity=sp.kronecker_product(z,z,z)
    require(even_a*fermion_parity==fermion_parity*even_a
            and even_b*fermion_parity==fermion_parity*even_b
            and (even_a*even_b)[0,0]==sp.I, "even Majorana phase convention")
    for length in range(7):
        for word in product((0,1),repeat=length):
            even,model=sp.eye(8),sp.eye(4)
            for letter in word:
                even*=(even_a,even_b)[letter]
                model*=alphabet[letter]
            require(even[0,0]==model[0,0],"even Majorana realization changed a word")
    positive_mark=(sp.eye(4)+b)/2
    q=sp.Symbol("q",positive=True)
    positive_response=(vacuum.H*a*sp.diag(1,q,q**2,q**3)*positive_mark*vacuum)[0]
    require(b.eigenvals()=={-1:2,1:2} and positive_mark.eigenvals()=={0:2,1:2}
            and (vacuum.H*positive_mark*vacuum)[0]==sp.Rational(1,2)
            and positive_response==sp.I*q/2, "the positive lift must change both mark responses")
    return {"exact_words_per_realization":checked,"AB":"I","BA":"-I",
            "ABAB":str((a*b*a*b)[0,0]),"ordered_after_intermediate_cap":"0",
            "literal_early_positive_mark":{"original_scalar_nonnegative":True,
                "coherent_B_has_negative_eigenvalue":True,"positive_projection_vacuum_response":"1/2",
                "positive_projection_ordered_response":str(positive_response)}}


def pair(value):
    return [str(sp.re(value)),str(sp.im(value))]


def form(matrix):
    return [[pair(value) for value in row] for row in matrix.tolist()]


def contract_candidate(q, empty_step, reference_step, orientation):
    """Compute the finite candidate's four marks and reflected blocks in the R113 reader."""
    a,b,parity=marked_operators()
    b=orientation*b
    t=sp.diag(*(empty_step*q**j for j in range(4)))
    cap=sp.diag(1,0,0,0)
    bare=t**3
    paired_word=a*b*a*b
    transfer=-sp.I*a*b
    reference=reference_step**3

    def leg(operator):
        return {"ambient_operator":form(operator)}

    def history(operator):
        return {"time_id":"oscillator-control","durations":[1,1,1],
                "early":leg(operator),"middle":leg(operator),"late":leg(operator)}

    with TemporaryDirectory(prefix="originlab-coherent-") as tmp:
        root=Path(tmp)
        origin_file=Path(__file__).resolve()
        refs=[{"path":str(origin_file),"sha256":digest(origin_file)}]
        domain={"domain_id":"coherent-oscillator-arithmetic","source_originals":refs}
        base=root/"model.json"
        base.write_text(json.dumps(domain))
        binding={"path":base.name,"sha256":digest(base)}
        bundle={
            "schema":"tfpt.original-blocked-source.v1","input_role":"arithmetic_control_only",
            "matrix_convention":"sesquilinear_forms_in_original_word_basis",
            "ensemble_scope":"whole_original_ensemble","domain_id":domain["domain_id"],
            "native_inputs_sha256":binding["sha256"],"source_refs":refs,
            "component_bindings":dict.fromkeys(COMPONENTS,"coherent oscillator control"),
            "retained_words":[["number",j] for j in range(4)],
            "time":{"kind":"original_euclidean_ordered","id":"oscillator-control","T":3,"t_U":1,"t_Y":2},
            "gram_form":form(sp.eye(4)),"physical_selector_form":form(cap),
            "primitive_ambient_basis":{"source_binding":"four oscillator number states",
                                       "gram_form":form(sp.eye(4)),"word_coordinates":form(sp.eye(4))},
            "boundary_pairing_form":form(cap),"scalar_factor":pair(1),
            "reference_window_history":{
                **history(reference_step*sp.eye(4)),
                "reference_binding":"independent scalar transfer on the vacuum reference",
                "selector_form":form(cap),"boundary_pairing_form":form(cap),"scalar_factor":pair(1),
            },
            "representation":"original_exact_finite_fiber_sum",
            "configuration_unit":{"native_map_binding":"single four state control fibre",
                                  "kind":"whole_original_history_with_both_bosonic_endpoints"},
            "measure_rule":{"kind":"original_finite_atomic_measure","source_binding":"unit point mass"},
            "fiber_records":[{
                "original_configuration_id":"four-state-control","whole_fiber_binding":"number states 0 through 3",
                "original_measure_binding":"unit point mass","original_measure_mass":pair(1),
                "unnormalized_scalar_weight":pair(1),
                "ordered_window_legs":{**history(t),"Y":leg(a),"U_operator":leg(b),"pair_transfer":leg(transfer)},
            }],
            "reflected_pair_cap":{
                "source_id":f"{domain['domain_id']}:{binding['sha256']}","time_id":"oscillator-control",
                "ordered_word":["A","B","A","B"],"paired_support_binding":"whole even/odd pair algebra",
                "operator_dictionary_binding":"oscillator A and B","two_window_order_binding":"control",
                "absolute_reference_binding":"independently chosen reference transfer",
                "paired_projector_form":form(sp.eye(4)),"ordered_abab_form":form(bare*paired_word),
                "same_history_defect":{
                    "original_mixed_moment_binding":"same control transfer and input",
                    "from_fiber_transfers":True,
                },
            },
        }
        attached=root/"run.json"
        attached.write_text(json.dumps({**domain,"original_domain_binding":binding,"whole_source_bundle":bundle}))
        result=evaluate(root,attached)
    z=sp.simplify((empty_step/reference_step)**3)
    expected={"Z0":str(z),"ZU":"0","ZY":"0","ZUY":str(sp.I*orientation*q*z)}
    require(result["values"]==expected and not result["physical_source_identity_proved"],
            "coherent receipt differs from direct ordered contraction")
    require(result["ordered_primitive_atom_count"]==1 and result["absolute_reference_from_primitive_history"],
            "coherent search bypassed primitive source or reference construction")
    probe=result["reflected_pair_cap"]
    require(probe["status"]=="EXACT_PAIRED_WORD_MATCH" and probe["normalized_ABAB"]=="-1",
            "the coherent model failed the ABAB control")
    defects=probe["same_history_defect"]
    require(defects["derived_from_fiber_transfers"] and defects["form_rank"]==2
            and defects["selected_rank"]==int(orientation<0),
            "quadratic Gram failed to distinguish the reversed phase after the cap")
    return {"q":str(q),"empty_step":str(empty_step),"reference_step":str(reference_step),
            "orientation":orientation,"reference":str(reference),"marks":result["values"],
            "connected_UY":result["connected_UY"],"ABAB":probe["normalized_ABAB"],
            "same_history_defect":defects}


def search(stage=2):
    """Enumerate one finite shell in the sum of three reduced denominators."""
    if type(stage) is not int or stage<2:
        raise ValueError("coherent stage must be an integer >= 2")
    rates={sp.Rational(n,d) for d in range(2,stage+1) for n in range(1,d)}
    # The normalized empty/reference factors are independent. Vary each
    # without fitting any physical target or claiming an absolute scale.
    vacua={sp.S.One,*rates}
    candidates=[(q,v,r) for q,v,r in product(sorted(rates),sorted(vacua),sorted(vacua))
                if q.q+v.q+r.q==stage+2]
    receipts=[contract_candidate(q,empty,reference,phase)
              for q,empty,reference in candidates for phase in (1,-1)]
    nominal=contract_candidate(sp.Rational(1,2),sp.S.One,sp.S.One,1)
    rate_alias=contract_candidate(sp.Rational(1,3),sp.S.One,sp.S.One,1)
    vacuum_alias=contract_candidate(sp.Rational(1,2),sp.Rational(1,2),sp.S.One,1)
    require(nominal["ABAB"]==rate_alias["ABAB"]==vacuum_alias["ABAB"]
            and nominal["connected_UY"]!=rate_alias["connected_UY"]
            and nominal["connected_UY"]==vacuum_alias["connected_UY"]
            and nominal["marks"]["Z0"]!=vacuum_alias["marks"]["Z0"],
            "clock and vacuum freedom require distinguishable exact source receipts")
    return {"status":"EXACT_COHERENT_CANDIDATES", "stage":stage,
            "candidate_count":len(receipts),"candidates":receipts,
            "gaussian":gaussian_covariance(),"number_composition":number_sectors(stage),
            "decoder":decoder_controls(),"ordered_words":words_control(),
            "non_identifiability":{"nominal":nominal,"changed_rate":rate_alias,
                                   "changed_empty_transfer":vacuum_alias},
            "input_role":"arithmetic_control_only",
            "physical_solution":False}
