"""Original R113 reader and exact contraction for explicitly supplied forms."""
from pathlib import Path
import json
from hashlib import sha256
import sympy as sp

from language import positive
from verify import digest

NATIVE = Path("experiments/theory-contracts/original-paper-source-reconciliation-20261102/"
              "contribution_20261106Q/domain/native-diagnostic-v1.json")
MARKS = ("Z0", "ZU", "ZY", "ZUY")
COMPONENTS = (
    "all_bosonic_fibres", "whole_physical_selector", "boson_generator_or_history",
    "fermion_pairing_and_higher_vertices", "reference_and_eta", "absolute_rest",
    "same_original_time_and_source_marks",
)


def source_input_key(domain,root=None):
    if not domain.is_file():
        return sha256(json.dumps(["missing_original_domain",str(domain)]).encode()).hexdigest()
    raw = domain.read_bytes()
    data = json.loads(raw)
    paths = []
    bundle = data.get("whole_source_bundle")
    if isinstance(bundle,str):
        packet_path = (domain.parent/bundle).resolve()
        paths.append(packet_path)
        if packet_path.is_file():
            try:
                bundle = json.loads(packet_path.read_bytes())
            except (ValueError,UnicodeError):
                bundle = None  # The reader will report the malformed packet.
    if binding := data.get("original_domain_binding"):
        paths.append((domain.parent/binding["path"]).resolve())
    if root is not None:
        references = list(data.get("source_originals",[]))
        if isinstance(bundle,dict) and isinstance(bundle.get("source_refs"),list):
            references.extend(bundle["source_refs"])
        paths.extend((root/ref["path"]).resolve() for ref in references
                     if isinstance(ref,dict) and isinstance(ref.get("path"),str))
    components = [sha256(raw).hexdigest(),*[digest(path) if path.is_file() else None for path in paths]]
    return sha256(json.dumps(components,separators=(",",":")).encode()).hexdigest()


def missing_source_result(root, error):
    if error.filename is None:
        return None
    root = root.resolve()
    path = Path(error.filename).resolve()
    if not any(path.is_relative_to(root/area) for area in ("docs/research","experiments/theory-contracts","verification")):
        return None
    return {"status":"BLOCKED","reason":"ORIGINAL_SOURCE_FILE_MISSING",
            "missing":[str(path.relative_to(root))],"physical_source_identity_proved":False}


def rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, str, sp.Rational)):
        raise ValueError("exact rational string or integer required")
    return sp.Rational(value)


def complex_pair(value):
    if not isinstance(value, list) or len(value) != 2:
        raise ValueError("scalar must be an exact [real, imaginary] pair")
    return rational(value[0])+sp.I*rational(value[1])


def as_form(data, size):
    if len(data) != size or any(len(row) != size for row in data):
        raise ValueError("form has the wrong word-basis dimensions")
    return sp.Matrix([[complex_pair(value) for value in row] for row in data])


def source_refs(root, records):
    if not records:
        raise ValueError("original input references are required")
    checked = []
    for record in records:
        path = Path(record["path"])
        path = (root/path).resolve()
        if digest(path) != record["sha256"]:
            raise ValueError(f"source hash mismatch for {path}")
        checked.append({"path":str(path.relative_to(root) if path.is_relative_to(root) else path),
                        "sha256":record["sha256"], "id":record.get("id")})
    return checked


def missing_reflected_pair_cap():
    return {"status":"BLOCKED", "reason":"REFLECTED_PAIR_CAP_PROBE_MISSING",
            "measured_ordered_numerator":None, "predicted_ordered_numerator":None,
            "normalized_ABAB":None, "paired_ratio":None,
            "missing":["bound paired support and independently measured ordered ABAB form on R113"]}


def blocked(domain, path, bindings):
    return dict(status="BLOCKED", domain=domain["domain_id"], domain_sha256=digest(path),
                source_refs=bindings, original_values=dict.fromkeys(MARKS), connected_UY=None,
                first_operation="evaluate_original_full_blocked_endpoint_source",
                missing=["original_blocked_fiber_packet", "physical whole weighted endpoint form",
                         "bound absolute empty reference"],
                fibre_pushforward={"status":"BLOCKED",
                    "missing":["two source-bound levels with weighted forms on every fibre",
                               "the common time, configuration map, word embedding and two absolute references"]},
                relational_original={"status":"BLOCKED","values":[],
                    "missing":["two physically admissible charged endpoint fields",
                               "two original regional transports with a common endpoint and charge convention",
                               "joint positive weighted forms and bound reference for those paths"]},
                four_mark_original={"status":"BLOCKED","values":[],
                    "missing":["four original oriented operations with squares and adjoints on one weighted source",
                               "original Cartan and charged loop dictionary on the same domain"]},
                reflected_pair_cap=missing_reflected_pair_cap(),
                ordered_return={"status":"BLOCKED","values":[],
                    "missing":["bound original charged insertion and adjoint on the reflected word space",
                               "bound positive time step and outer physical selector on that same space"]},
                energy_moments={"status":"BLOCKED",
                    "missing":["first and second moments of the same original generator on a common reflected word space"]},
                affine_polarization={"status":"BLOCKED",
                    "missing":["weighted bare and Cartan-prepared three-leg words at shared ordered times",
                               "same original preparation, adjoint and absolute empty reference bindings"]},
                affine_profile_returns={"status":"BLOCKED",
                    "missing":["weighted neutral grade-zero and full-root charged grade-one words",
                               "same original profile product, adjoint and absolute empty reference bindings"]},
                weak_split_returns={"status":"BLOCKED",
                    "missing":["both ordered weak-odd splits of the marked middle root on the weighted word space"]},
                weak_gauge_responses={"status":"BLOCKED",
                    "missing":["native local gauge dressing and reflected words for the chosen gauge singlet triangles",
                               "early and late transport and Casimir insertions on the same weighted form",
                               "early weak and quarter coweight insertion for fixed-middle path reconstruction"]},
                weak_clock_sectors={"status":"BLOCKED",
                    "missing":["the same early joint action at powers zero through seven on one original weighted endpoint and time dictionary"]},
                p1_charged_jet={"status":"BLOCKED",
                    "missing":["the whole two-face covariance, anchored hypercharge and heat jets on the original source"]},
                weak_parity_bridge={"status":"BLOCKED",
                    "missing":["original alpha Cartan, fermion parity and charged endpoint field on a joint reflected basis"]},
                weak_return_grade_two={"status":"BLOCKED",
                    "missing":["weighted grade-two excess words for the 18 fixed-middle weak legs on one reflected source"]},
                obligations={
                    "original_configuration_map":"Produce the ordered original U/Y/UY forms at both endpoints.",
                    "weight":"Derive the full original boson measure, fermion weight, eta and rest on the same domain.",
                    "selection":"Contract with the single original outer admissibility projector and cap.",
                    "reference":"Derive the bound empty reference with the same regulator and time dictionary.",
                    "limit":"Establish the joint original regulator and all-field limits after finite contraction."
                })


def fibre_forms(atom, size, u_time, y_time):
    if "ordered_window_kernels" not in atom:
        return {k:as_form(atom["ordered_source_forms"][k],size) for k in MARKS}
    kernel = atom["ordered_window_kernels"]
    if ("ordered_source_forms" in atom or not u_time<y_time or
            kernel.get("kind")!="whole_history_fixed_bosonic_scalar_U" or
            not isinstance(kernel.get("original_scalar_mark_binding"),str) or
            not kernel["original_scalar_mark_binding"].strip()):
        raise ValueError("two kernels require one bound scalar early mark on this strictly ordered history")
    u = rational(kernel["U"])
    if u<0:
        raise ValueError("original squared Higgs mark must be nonnegative")
    k0,ky = (as_form(kernel[k],size) for k in ("K0","KY"))
    return dict(Z0=k0,ZU=u*k0,ZY=ky,ZUY=u*ky)


def iter_weighted_fibres(atoms, size, u_time, y_time):
    """Integrate supplied oriented atom weights without discarding their identities."""
    seen = set()
    for atom in atoms:
        identifier = atom["original_configuration_id"]
        if not identifier or identifier in seen or not atom["whole_fiber_binding"] or not atom["original_measure_binding"]:
            raise ValueError("fibre identities and original bindings must be distinct and complete")
        seen.add(identifier)
        mass = complex_pair(atom["original_measure_mass"])
        if not mass.is_nonnegative:
            raise ValueError("atomic measure mass must be nonnegative real")
        weight = mass*complex_pair(atom["unnormalized_scalar_weight"])
        yield identifier, {mark:weight*form for mark,form in fibre_forms(atom,size,u_time,y_time).items()}


def fibre_pushforward_data(data, fine_labels, fine_size, time_id, u_time, y_time, fine_scale):
    """Read the common map, transported coarse forms and absolute scale."""
    rows = data["word_embedding"]
    coarse_size = len(rows[0]) if rows else 0
    if (not coarse_size or len(rows) != fine_size or any(len(row) != coarse_size for row in rows)
            or not data["configuration_map_binding"] or not data["word_embedding_binding"]
            or not data["reference_binding"] or data["time_id"] != time_id):
        raise ValueError("the fibre comparison requires the bound configuration map, word embedding, time and reference")
    embedding = sp.Matrix([[complex_pair(entry) for entry in row] for row in rows])
    if embedding.rank() != coarse_size:
        raise ValueError("the coarse word embedding must be injective")
    coarse = dict(iter_weighted_fibres(data["coarse_fiber_records"],coarse_size,u_time,y_time))
    assignment = data["fine_to_coarse"]
    if (not coarse or not isinstance(assignment,dict) or set(assignment) != set(fine_labels)
            or any(not isinstance(name,str) for name in assignment.values())
            or set(assignment.values()) != set(coarse)):
        raise ValueError("the fine configuration map must cover both supplied levels")
    reference = complex_pair(data["coarse_reference_scalar"])
    scalar = complex_pair(data["coarse_scalar_factor"])
    if not reference or not scalar or not fine_scale:
        raise ValueError("both levels need a nonzero absolute reference and weight scale")
    return embedding,coarse,assignment,scalar/reference


def check_fibre_pushforward(data, fine, fine_size, time_id, u_time, y_time, fine_scale):
    """Compare the supplied marked forms before the physical endpoint contraction."""
    embedding,coarse,assignment,coarse_scale = fibre_pushforward_data(
        data,fine,fine_size,time_id,u_time,y_time,fine_scale)
    coarse_size = embedding.cols
    pushed = {name:{k:sp.zeros(coarse_size) for k in MARKS} for name in coarse}
    for name, forms in fine.items():
        for mark, form in forms.items():
            pushed[assignment[name]][mark] += embedding.H*form*embedding
    for name, forms in coarse.items():
        for mark, form in forms.items():
            delta = (fine_scale*pushed[name][mark] - coarse_scale*form).applyfunc(sp.simplify)
            if not delta.is_zero_matrix:
                return {"status":"REFUTED_SUPPLIED_MARKED_PUSHFORWARD","first_coarse_configuration":name,
                        "first_mark":mark,"defect_rank":delta.rank()}
    return {"status":"EXACT_SUPPLIED_MARKED_PUSHFORWARD",
            "fine_fibres":len(fine),"coarse_fibres":len(coarse),"checked_forms":len(coarse)*len(MARKS)}


def path_contraction(data, marks, time, source_id, size, contract):
    from reconstruction import require
    from relational import network, path_angles, phase, zero
    require(data.get("source_id")==source_id and data.get("time_id")==time["id"],
            "regional forms must bind this same native source and time dictionary")
    operations = data["paths"]
    require(isinstance(operations,list) and len(operations)>=2,
            "at least two ordered regional path forms required")
    charge = data["endpoint_charge"]
    require(type(charge) is int and charge!=0 and isinstance(data.get("charge_spin_binding"),str) and
            data["charge_spin_binding"].strip(),"original charge and spin dictionary required")
    for key in ("endpoint_field_binding","transport_binding","admissibility_binding"):
        require(isinstance(data.get(key),str) and data[key].strip(),f"missing {key}")
    comparison = network(data["connection"])
    path_pairs = comparison["path_comparisons"]
    require(len(operations)==2 and len(path_pairs)==1,"use exactly two paths for a same-endpoint comparison")
    names = [item["id"] for item in operations]
    require(all(isinstance(name,str) and name for name in names) and len(set(names))==2,
            "distinct named regional paths required")
    for item,route in zip(operations,data["connection"]["path_pairs"][0]["routes"]):
        require(item.get("route")==route and isinstance(item["word"],list) and item["word"],
                "a bound regional word and ordered connection path are required for each form")
    require(operations[0]["word"]!=operations[1]["word"],"original path words must be distinct")
    path_forms = [as_form(item["form"],size) for item in operations]
    if "reflected_gram" in data:
        reflected = as_form(data["reflected_gram"],2)
        require(positive(reflected),"supplied reflected two-path Gram must be positive")
    first,second = (contract(form) for form in path_forms)
    empty = marks["Z0"]
    result = {"status":"EXACT_SUPPLIED_REGIONAL_CONTRACTION","source_id":source_id,"time_id":time["id"],
              "paths":names,"endpoint_charge_in_supplied_dictionary":charge,
              "absolute_responses":[str(first),str(second)],
              "relative_to_same_empty":[str(first/empty),str(second/empty)] if empty else None,
              "response_ratio":str(sp.simplify(first/second)) if second else None,
              "path_transport":path_pairs[0],
              "supplied_reflected_path_gram_psd":"reflected_gram" in data,
              "physical_original_fields_or_connection_derived":False}
    if "unit_charge_path_comparison_binding" in data:
        binding = data["unit_charge_path_comparison_binding"]
        require(isinstance(binding,str) and binding.strip() and data["connection"].get("dimension",1)==1,
                "Abelian response comparison needs a bound unit-charge convention")
        paths = data["connection"]["path_pairs"][0]["routes"]
        angle_a,angle_b = path_angles(data["connection"]["edges"],paths)
        expected = phase(charge*(angle_a-angle_b))
        residual = sp.simplify(first-expected*second)
        result["response_holonomy_check"] = {
            "binding":binding,"expected_character":str(expected),"residual":str(residual),
            "status":"DISTINGUISHED" if not zero(residual) else "EXACT_MATCH" if second else "ZERO_RESPONSES"}
    return result,dict(zip(names,(first,second)))


def joint_source(data):
    """Exact connected responses and declared linear Ward identities."""
    from reconstruction import scalar
    marks = data["responses"]
    return check_joint_marks(data,{name:scalar(value) for name,value in marks.items()})


def check_joint_marks(data,values):
    from reconstruction import require, scalar
    require(values.get("Z0") not in (None,0),"common empty response must be nonzero")
    response = []
    for row in data["mixed"]:
        a,b,ab = (row[k] for k in ("first","second","joint"))
        value = values[ab]/values["Z0"]-values[a]*values[b]/values["Z0"]**2
        response.append({"first":a,"second":b,"joint":ab,"connected":str(value)})
    ward = {}
    for row in data.get("ward_identities",[]):
        name = row["id"]
        require(isinstance(name,str) and name and name not in ward,"distinct Ward identities required")
        defect = sum(scalar(c)*values[k] for k,c in row["coefficients"].items())-scalar(row["rhs"])
        ward[name] = {"satisfied":defect==0,"exact_residual":str(defect)}
    return {"status":"EXACT_SUPPLIED_COMMON_MARKS","mixed":response,"linear_identities":ward,
            "scope":"supplied common ordered responses and supplied linear identities"}


def evaluate(root, path=None):
    path = path or root/NATIVE
    domain = json.loads(path.read_text())
    bindings = source_refs(root, domain["source_originals"])
    bundle = domain.get("whole_source_bundle")
    if bundle is None:
        return blocked(domain, path, bindings)
    packet_hash = None
    if isinstance(bundle, str):
        packet_path = (path.parent/bundle).resolve()
        if not packet_path.is_file():
            result = blocked(domain,path,bindings)
            result["missing_external_packet"] = str(packet_path)
            return result
        packet_hash = digest(packet_path)
        bundle = json.loads(packet_path.read_text())
    role = bundle["input_role"]
    if role not in ("source_reviewed_original", "arithmetic_control_only"):
        raise ValueError("unknown bundle role")
    required = dict(schema="tfpt.original-blocked-source.v1",
                    matrix_convention="sesquilinear_forms_in_original_word_basis",
                    ensemble_scope="whole_original_ensemble",
                    domain_id=domain["domain_id"])
    if any(bundle.get(key) != value for key,value in required.items()):
        raise ValueError("incompatible original packet convention or domain")
    # External bundles bind to the original domain before the bundle is attached.
    # For inline or attached packets the domain names that immutable input explicitly.
    binding = domain.get("original_domain_binding")
    if not isinstance(binding, dict):
        raise ValueError("supply original_domain_binding with original native input path and SHA256")
    original = source_refs(path.parent, [binding])[0]
    original_path = (path.parent/original["path"]).resolve()
    original_domain = json.loads(original_path.read_text())
    if original_domain.get("domain_id") != domain["domain_id"]:
        raise ValueError("original native domain ID mismatch")
    domain_inputs = {k:v for k,v in domain.items() if k not in {"original_domain_binding","whole_source_bundle"}}
    native_inputs = {k:v for k,v in original_domain.items() if k != "whole_source_bundle"}
    if domain_inputs != native_inputs or bundle["native_inputs_sha256"] != binding["sha256"]:
        raise ValueError("packet or declared inputs do not bind the same native domain")
    if any(not isinstance(bundle["component_bindings"].get(key),str)
           or not bundle["component_bindings"][key].strip() for key in COMPONENTS):
        raise ValueError("one or more physical component bindings are missing")
    refs = source_refs(root, bundle["source_refs"])
    words = bundle["retained_words"]
    size = len(words)
    if not size or any(not isinstance(w,list) for w in words) or len({json.dumps(w) for w in words}) != size:
        raise ValueError("original word basis must contain distinct explicit labels")
    time = bundle["time"]
    total, u_time, y_time = (rational(time[k]) for k in ("T","t_U","t_Y"))
    if time["kind"] != "original_euclidean_ordered" or not time["id"] or not 0 <= u_time <= y_time <= total:
        raise ValueError("original Euclidean time ordering is missing or invalid")
    if domain.get("source_time_id",time["id"]) != time["id"]:
        raise ValueError("packet and domain use different time dictionaries")
    gram = as_form(bundle["gram_form"],size)
    if not positive(gram):
        raise ValueError("original word Gram is not Hermitian PSD")
    _, pivots = gram.rref()
    if not pivots:
        raise ValueError("empty OS quotient")
    null = gram.nullspace()
    quotient_gram = gram.extract(pivots,pivots)
    metric = quotient_gram.inv()

    def quotient(form):
        if any(form*v != sp.zeros(size,1) or v.H*form != sp.zeros(1,size) for v in null):
            raise ValueError("a supplied form does not descend to the original OS quotient")
        return metric*form.extract(pivots,pivots)

    def projector(value):
        form = as_form(value,size)
        operator = quotient(form)
        if form != form.H or operator*operator != operator:
            raise ValueError("selector is not an orthogonal projector in the original metric")
        return operator

    ambient = bundle.get("primitive_ambient_basis")
    if ambient is not None:
        if not isinstance(ambient.get("source_binding"),str) or not ambient["source_binding"].strip():
            raise ValueError("ambient basis needs an original coordinate binding")
        ambient_size = len(ambient["gram_form"])
        ambient_gram = as_form(ambient["gram_form"],ambient_size)
        data = ambient["word_coordinates"]
        if len(data) != ambient_size or any(len(row) != size for row in data):
            raise ValueError("ambient word coordinates have the wrong dimensions")
        coordinates = sp.Matrix([[complex_pair(value) for value in row] for row in data])
        if not ambient_size or not positive(ambient_gram) or coordinates.H*ambient_gram*coordinates != gram:
            raise ValueError("ambient coordinates do not reproduce the retained reflected Gram")
        ambient_null = ambient_gram.nullspace()

    def ambient_operator(data):
        if ambient is None or "form" in data or "image_norm_form" in data:
            raise ValueError("ambient operator needs one bound coordinate basis and one representation")
        operator = as_form(data["ambient_operator"],ambient_size)
        if any(ambient_gram*operator*word != sp.zeros(ambient_size,1) for word in ambient_null):
            raise ValueError("primitive operator does not descend to the ambient reflected quotient")
        return operator

    def ambient_form(operator):
        return coordinates.H*ambient_gram*operator*coordinates

    def leg_forms(data):
        if "ambient_operator" in data:
            image = ambient_operator(data)*coordinates
            form = coordinates.H*ambient_gram*image
            image_norm = image.H*ambient_gram*image
        else:
            if not isinstance(data.get("image_norm_binding"),str) or not data["image_norm_binding"].strip():
                raise ValueError("each primitive leg needs an independently sourced image norm")
            form = as_form(data["form"],size)
            image_norm = as_form(data["image_norm_form"],size)
        return form,image_norm

    def closed_leg(data):
        form,image_norm = leg_forms(data)
        operator = quotient(form)
        quotient(image_norm)
        if operator.H*quotient_gram*operator != image_norm.extract(pivots,pivots):
            raise ValueError("primitive leg does not close on the retained reflected word space")
        return operator

    def window_legs(history,reader):
        if (history.get("time_id") != time["id"]
                or tuple(map(rational,history["durations"])) != (u_time,y_time-u_time,total-y_time)):
            raise ValueError("primitive history and original window times differ")
        return tuple(reader(history[key]) for key in ("early","middle","late"))

    def marked_history(history,lift):
        if not u_time < y_time or ("U" in history) == ("U_operator" in history):
            raise ValueError("specify one U insertion and strictly ordered marked windows")
        keys = ["early","middle","late","Y"]
        if "U_operator" in history:
            keys.append("U_operator")
        use_ambient = all("ambient_operator" in history[key] for key in keys)
        reader = ambient_operator if use_ambient else closed_leg
        early,middle,late = window_legs(history,reader)
        if "U_operator" in history:
            if role != "arithmetic_control_only":
                raise ValueError("original U must be the nonnegative scalar insertion of this fibre")
            u_early = reader(history["U_operator"])*early
        else:
            u = rational(history["U"])
            if u < 0:
                raise ValueError("original U insertion must be nonnegative")
            u_early = u*early
        late_y = late*reader(history["Y"])
        encode = ambient_form if use_ambient else lift
        return dict(Z0=encode(late*middle*early),ZU=encode(late*middle*u_early),
                    ZY=encode(late_y*middle*early),ZUY=encode(late_y*middle*u_early))

    selector = projector(bundle["physical_selector_form"])
    cap = quotient(as_form(bundle["boundary_pairing_form"],size))
    normalization = bundle.get("normalization_scope","absolute")
    if normalization not in ("absolute","relative_mark_contraction","relative_radial_moment_bound"):
        raise ValueError("unknown original normalization scope")
    scalar = complex_pair(bundle["scalar_factor"])
    radial_only = normalization=="relative_radial_moment_bound"
    reference_history = bundle.get("reference_window_history")
    has_reference_scalar = bundle.get("reference_scalar") is not None
    if normalization!="absolute":
        if has_reference_scalar or reference_history is not None:
            raise ValueError("use the absolute path when an absolute reference is supplied")
        reference = None
    elif reference_history is not None:
        reference_binding = reference_history.get("reference_binding")
        if has_reference_scalar or not isinstance(reference_binding,str) or not reference_binding.strip():
            raise ValueError("supply one absolute reference with an independent empty history binding")
        use_ambient = all("ambient_operator" in reference_history[key] for key in ("early","middle","late"))
        early,middle,late = window_legs(reference_history,ambient_operator if use_ambient else closed_leg)
        empty_word = late*middle*early
        if use_ambient:
            empty_word = quotient(ambient_form(empty_word))
        ref_selector = projector(reference_history["selector_form"])
        ref_cap = quotient(as_form(reference_history["boundary_pairing_form"],size))
        ref_scalar = complex_pair(reference_history["scalar_factor"])
        reference = sp.simplify(ref_scalar*(ref_cap*ref_selector*empty_word*ref_selector).trace())
    else:
        reference = complex_pair(bundle["reference_scalar"])
    if reference == 0:
        raise ValueError("empty reference is zero")
    representation = bundle["representation"]
    pushforward = bundle.get("coarse_graining")
    if pushforward is not None and (representation != "original_exact_finite_fiber_sum" or reference is None):
        raise ValueError("the two-level comparison needs explicit fine fibre records")
    kernel_atoms = 0
    fibre_audit = None
    primitive_atoms = 0
    derived_defect = None
    if radial_only and (representation!="original_exact_blocked_endpoint_source"
                        or "radial_moment_probe" not in bundle
                        or set(bundle["ordered_source_forms"])!={"Z0","ZU"}):
        raise ValueError("radial-bound-only scope needs the two radial base forms and the bound moment probe")
    if representation == "original_exact_blocked_endpoint_source":
        form_keys = ("Z0","ZU") if radial_only else MARKS
        forms = {k:as_form(bundle["ordered_source_forms"][k],size) for k in form_keys}
        atoms = None
    elif representation == "original_exact_finite_fiber_sum":
        forms = {k:sp.zeros(size) for k in MARKS}
        atoms = bundle["fiber_records"]
        unit = bundle["configuration_unit"]
        measure = bundle["measure_rule"]
        if (not atoms or not unit.get("native_map_binding") or unit.get("kind") not in {
            "whole_original_history_with_both_bosonic_endpoints",
            "original_retained_field_fiber_with_whole_history_endpoint_packet",
        } or measure.get("kind") != "original_finite_atomic_measure" or not measure.get("source_binding")):
            raise ValueError("no original atomic domain, map or measure supplied")
        seen = set()
        fibres = {}
        basis_coordinates = metric*gram.extract(pivots,range(size))

        def lift(operator):
            return basis_coordinates.H*quotient_gram*operator*basis_coordinates

        defect_data = bundle.get("reflected_pair_cap",{}).get("same_history_defect",{})
        if defect_data.get("from_fiber_transfers") is True:
            if "form" in defect_data or "moments" in defect_data or reference is None:
                raise ValueError("derived same-history defect requires its own input and an absolute reference")
            derived_defect = [sp.zeros(size) for _ in range(3)]
        for atom in atoms:
            identifier = atom["original_configuration_id"]
            if not identifier or identifier in seen or not atom["whole_fiber_binding"] or not atom["original_measure_binding"]:
                raise ValueError("fibre identities and original bindings must be distinct and complete")
            seen.add(identifier)
            mass = complex_pair(atom["original_measure_mass"])
            if not mass.is_nonnegative:
                raise ValueError("atomic measure mass must be nonnegative real")
            weight = mass*complex_pair(atom["unnormalized_scalar_weight"])
            history = atom.get("ordered_window_legs")
            if history is not None:
                if "ordered_source_forms" in atom or "ordered_window_kernels" in atom:
                    raise ValueError("primitive legs and contracted marks cannot both specify a fibre")
                atom_forms = marked_history(history,lift)
                primitive_atoms += 1
            else:
                atom_forms = fibre_forms(atom,size,u_time,y_time)
                kernel_atoms += int("ordered_window_kernels" in atom)
            if derived_defect is not None:
                if not weight.is_nonnegative or history is None or "pair_transfer" not in history:
                    raise ValueError("same-history Gram needs a positive weight and a transfer on every fibre")
                overlap,image_norm = leg_forms(history["pair_transfer"])
                delta = (gram,overlap,image_norm)
                for index,block in enumerate(delta):
                    derived_defect[index] += weight*block
            weighted = {mark:weight*form for mark,form in atom_forms.items()}
            for k,q in weighted.items():
                forms[k] += q
            if pushforward is not None:
                fibres[identifier] = weighted
        if pushforward is not None:
            fibre_audit = check_fibre_pushforward(pushforward,fibres,size,time["id"],u_time,y_time,scalar/reference)
            if fibre_audit["status"] == "REFUTED_SUPPLIED_MARKED_PUSHFORWARD":
                return {"status":fibre_audit["status"],"input_role":role,
                        "domain":domain["domain_id"],"source_refs":bindings+refs,
                        "fibre_pushforward":fibre_audit,"original_values":dict.fromkeys(MARKS),
                        "physical_source_identity_proved":False,"physical_limits_proved":False}
    else:
        raise ValueError("unsupported endpoint representation")
    def operator_numerator(op):
        return sp.simplify(scalar*(cap*selector*op*selector).trace())

    def numerator(form):
        return operator_numerator(quotient(form))

    raw = {k:numerator(q) for k,q in forms.items()}
    marks = {k:sp.simplify(value/reference) for k,value in raw.items()} if reference is not None else {}
    z = raw["Z0"]
    connected = sp.simplify(raw["ZUY"]/z-raw["ZU"]*raw["ZY"]/z**2) if z and not radial_only else None
    source_id = f"{domain['domain_id']}:{binding['sha256']}"
    candidate = dict(status="EXACT_SUPPLIED_CONTRACTION", input_role=role, domain=domain["domain_id"],
                domain_sha256=digest(path), native_domain_sha256=binding["sha256"],
                packet_sha256=packet_hash, source_refs=bindings+refs, quotient_rank=len(pivots),
                time=time, atom_count=len(atoms) if atoms is not None else None,
                ordered_two_kernel_atom_count=kernel_atoms,
                ordered_primitive_atom_count=primitive_atoms,
                absolute_reference_from_primitive_history=reference_history is not None,
                absolute_reference=str(reference) if reference is not None else None,
                unnormalized_numerators={k:str(v) for k,v in raw.items()},
                values={k:str(v) for k,v in marks.items()},
                original_values={k:str(marks[k]) for k in MARKS} if marks and role=="source_reviewed_original"
                else dict.fromkeys(MARKS),
                connected_UY=str(connected) if connected is not None else None,
                physical_source_identity_proved=False, physical_limits_proved=False)
    if "radial_moment_probe" in bundle:
        from r113_reflected_radial import RADIAL_DATA, RADIAL_PROFILE_SHA256, moment_bound

        probe = bundle["radial_moment_probe"]
        if (not isinstance(probe,dict) or probe.get("source_id")!=source_id
                or probe.get("time_id")!=time["id"] or not z.is_positive
                or probe.get("radial_profile_sha256")!=RADIAL_PROFILE_SHA256
                or digest(root/RADIAL_DATA)!=RADIAL_PROFILE_SHA256
                or any(not isinstance(probe.get(key),str) or not probe[key].strip()
                       for key in ("positive_radial_law_binding","normalized_late_profile_binding",
                                   "paired_radial_operator_binding"))):
            raise ValueError("radial moment probe needs the same positive source, time and reviewed profile")
        probe_forms = ("late_radial_form","early_radial_square_form","late_radial_square_form","cross_radial_form")
        relative_moments = [sp.simplify(numerator(as_form(probe[key],size))/z) for key in probe_forms]
        candidate["radial_moment_probe"] = {
            **moment_bound(sp.simplify(raw["ZU"]/z),*relative_moments,connected=connected),
            "source_id":source_id,"time_id":time["id"],"normalization":"common_empty",
            "scope":"conditional on the supplied radial law and normalized Y profile bindings"}
    if reference is None:
        status = "EXACT_RADIAL_MOMENT_BOUND" if radial_only else "EXACT_RELATIVE_CONTRACTION"
        missing = ["bound absolute empty reference"]
        if radial_only:
            missing.append("original numeric Y and UY marks")
        return {**candidate,"status":status if z else "ZERO_EMPTY_CONTRACTION",
                "original_values":dict.fromkeys(MARKS),"absolute_reference":None,
                "normalized_by_common_empty":{k:str(sp.simplify(v/z)) for k,v in raw.items()} if z else {},
                "relative_observables_only":True,"missing":missing}
    if "q0_rule_probe" in bundle:
        probe = bundle["q0_rule_probe"]
        if (probe.get("source_id") != source_id or probe.get("time_id") != time["id"]
                or any(not isinstance(probe.get(key), str) or not probe[key].strip()
                       for key in ("operator_dictionary_binding", "two_window_time_binding",
                                   "absolute_reference_binding"))
                or not total > 0 or 2*rational(probe["step"]) != total):
            raise ValueError("Q0 probe requires the same original source, two windows and absolute reference")
        endpoints = []
        for side in ("left", "right"):
            index = probe[f"{side}_index"]
            word = probe[f"{side}_current_word"]
            if (type(index) is not int or not 0 <= index < size or not isinstance(word, list)
                    or any(not isinstance(letter, list) or len(letter) != 2
                           or letter[0] not in ("H", "E", "F")
                           or type(letter[1]) is not int or letter[1] >= 0 for letter in word)):
                raise ValueError("invalid Q0 affine word endpoint")
            vector = metric*gram.extract(pivots, [index])
            if selector*vector != vector or gram[index, index] == 0:
                raise ValueError("Q0 endpoint must survive the original outer selector")
            endpoints.append((vector, gram[index, index]))
        (left, left_norm), (right, right_norm) = endpoints
        reduced_gram = gram.extract(pivots, pivots)
        value = sp.simplify(scalar*(left.H*reduced_gram*selector*quotient(forms["Z0"])*selector*right)[0])
        unit = sp.sqrt(left_norm*right_norm)
        candidate["q0_rule_probe"] = {
            "input_role": role, "source_id": source_id, "time_id": time["id"],
            "left": probe["left_current_word"], "right": probe["right_current_word"],
            "absolute_numerator": str(value), "absolute_reference": str(reference),
            "unit_word_response": str(sp.simplify(value/reference/unit)),
        }
    if "composed_search" in bundle:
        from composed_search import run
        data = bundle["composed_search"]
        if (data.get("source_id") != source_id or data.get("time_id") != time["id"] or
                data.get("reference_binding") != bundle["component_bindings"]["reference_and_eta"]):
            raise ValueError("composed search must bind the original source, time and reference")
        labels = data["mark_words"]
        if set(labels) != set(MARKS):
            raise ValueError("composed search must account for all four original marks")
        def observation(word, value):
            real, imaginary = value.as_real_imag()
            return {"word":word,"value":[str(real),str(imaginary)]}
        observations = {"fit":[observation(labels[k],marks[k]) for k in MARKS if k != "ZUY"],
                        "withheld":[observation(labels["ZUY"],marks["ZUY"])]}
        for row in data.get("withheld_forms", []):
            observations["withheld"].append(observation(row["word"],numerator(as_form(row["form"],size))/reference))
        candidate["composed_search"] = {
            **run(data["case"],gram.extract(pivots,pivots),(scalar/reference)*selector*cap*selector,observations),
            "source_id":source_id,"time_id":time["id"],"input_role":role,
            "original_word_indices":list(pivots),
        }
    if fibre_audit is not None:
        candidate["fibre_pushforward"] = fibre_audit
    candidate["relational_original"] = {"status":"BLOCKED","values":[],
        "missing":["bound same-source regional endpoint fields, connection and two ordered path forms"]}
    candidate["four_mark_original"] = {"status":"BLOCKED","values":[],
        "missing":["bound original orientation, state, four ordered operations and integral charge dictionary"]}
    candidate["reflected_pair_cap"] = missing_reflected_pair_cap()
    if "reflected_pair_cap" in bundle:
        data = bundle["reflected_pair_cap"]
        if (not isinstance(data,dict) or data.get("source_id") != source_id
                or data.get("time_id") != time["id"] or data.get("ordered_word") != ["A","B","A","B"]
                or any(not isinstance(data.get(key),str) or not data[key].strip() for key in
                       ("paired_support_binding","operator_dictionary_binding","two_window_order_binding",
                        "absolute_reference_binding"))):
            raise ValueError("reflected pair cap requires the same original source, ordered ABAB, time and reference")
        paired_form = as_form(data["paired_projector_form"],size)
        paired = quotient(paired_form)
        if paired_form != paired_form.H or paired*paired != paired:
            raise ValueError("reflected pair cap requires a metric-orthogonal paired projector")
        measured = quotient(as_form(data["ordered_abab_form"],size))
        bare = quotient(forms["Z0"])
        commutator = paired*bare-bare*paired
        measured_numerator = operator_numerator(measured)

        def matrix_rows(op):
            return [[str(entry) for entry in row] for row in op.tolist()]

        restricted = {}
        for key in ("ZY","ZUY"):
            restricted_op = paired*quotient(forms[key])*paired
            restricted[key] = {"rank":restricted_op.rank(),
                               "selected_numerator":str(operator_numerator(restricted_op))}
        probe = {
            "status":"PAIR_CAP_DOES_NOT_REDUCE_BARE_TRANSFER",
            "input_role":role, "source_id":source_id, "time_id":time["id"],
            "quotient_basis_indices":list(pivots),
            "absolute_reference":str(reference), "bare_numerator":str(z),
            "paired_projector_rank":paired.rank(), "bare_commutator_rank":commutator.rank(),
            "measured_ordered_numerator":str(measured_numerator), "predicted_ordered_numerator":None,
            "measured_ordered_value_over_reference":str(sp.simplify(measured_numerator/reference)),
            "predicted_ordered_value_over_reference":None,
            "normalized_ABAB":str(sp.simplify(measured_numerator/z)) if z else None,
            "paired_ratio":None, "restricted_first_responses":restricted,
            "measured_operator":matrix_rows(measured), "predicted_operator":None,
            "matrix_residual":None, "selected_matrix_residual_rank":None,
            "physical_source_identity_proved":False,
        }
        if "same_history_defect" in data:
            defect = data["same_history_defect"]
            if (not isinstance(defect,dict)
                    or not isinstance(defect.get("original_mixed_moment_binding"),str)
                    or not defect["original_mixed_moment_binding"].strip()):
                raise ValueError("pairing defect requires the original same-history squared moment")
            if sum(("form" in defect,"moments" in defect,defect.get("from_fiber_transfers") is True)) != 1:
                raise ValueError("supply one squared defect form or its same-history Gram blocks")
            if defect.get("from_fiber_transfers") is True and derived_defect is None:
                raise ValueError("same-history fibre transfers require a finite atomic representation")
            if "moments" in defect or derived_defect is not None:
                blocks = derived_defect if derived_defect is not None else [
                    as_form(defect["moments"][name],size) for name in
                    ("input_norm","transfer_overlap","transfer_norm")]
                initial, overlap, transferred = blocks
                if not positive(initial.row_join(overlap).col_join(overlap.H.row_join(transferred))):
                    raise ValueError("same-history input and transfer blocks lack a common positive Gram")
                for form in blocks:
                    quotient(form)
                defect_form = transferred-overlap-overlap.H+initial
            else:
                defect_form = as_form(defect["form"],size)
            if not positive(defect_form):
                raise ValueError("same-history squared pairing defect must be Hermitian PSD")
            defect_op = quotient(defect_form)
            form_rank = defect_op.rank()
            selected_rank = (selector*defect_op*selector).rank()
            defect_numerator = operator_numerator(defect_op)
            probe["same_history_defect"] = {
                "derived_from_fiber_transfers":derived_defect is not None,
                "form_rank":form_rank,"selected_rank":selected_rank,
                "numerator":str(defect_numerator),
                "value_over_reference":str(sp.simplify(defect_numerator/reference)),
                "nonzero_hidden_by_selector":bool(form_rank and not selected_rank),
                "nonzero_selected_hidden_by_contraction":bool(selected_rank and not defect_numerator),
            }
        if commutator == sp.zeros(paired.rows):
            paired_op = paired*bare*paired
            paired_numerator = operator_numerator(paired_op)
            predicted = bare-2*paired_op
            residual = measured-predicted
            selected_residual = selector*residual*selector
            predicted_numerator = operator_numerator(predicted)
            probe.update(status="EXACT_PAIRED_WORD_MATCH" if residual == sp.zeros(paired.rows)
                         else "PAIRED_WORD_MATRIX_MISMATCH",
                         paired_numerator=str(paired_numerator),
                         paired_ratio=str(sp.simplify(paired_numerator/z)) if z else None,
                         predicted_ordered_numerator=str(predicted_numerator),
                         predicted_ordered_value_over_reference=str(sp.simplify(predicted_numerator/reference)),
                         contracted_residual=str(sp.simplify(measured_numerator-predicted_numerator)),
                         predicted_operator=matrix_rows(predicted),
                         matrix_residual=matrix_rows(residual),
                         selected_matrix_residual_rank=selected_residual.rank())
        candidate["reflected_pair_cap"] = probe
    if "four_mark_selection" in bundle:
        from selection import identify_face_refinement, operator_checks
        data = bundle["four_mark_selection"]
        if data.get("source_id")!=source_id or data.get("time_id")!=time["id"] or any(
            not isinstance(data.get(k),str) or not data[k].strip()
            for k in ("orientation_binding","adjoint_binding","state_binding")):
            raise ValueError("four operations require the same native source, time, and original orientation, adjoint and state")
        empty = data["identity_index"]
        if type(empty) is not int or not 0<=empty<size or words[empty]!=[]:
            raise ValueError("empty word must bind the same original state")
        labels = data["operations"]
        if not isinstance(labels,dict) or len(labels)!=4 or not all(isinstance(k,str) and k for k in labels):
            raise ValueError("four named original oriented operations required")
        ops = [quotient(as_form(q,size)) for q in labels.values()]
        reduced_metric = gram.extract(pivots,pivots)
        omega = metric*gram.extract(pivots,[empty])
        transport_inputs = charge_inputs = None
        if "transport" in data:
            transport = data["transport"]
            if not isinstance(transport.get("binding"),str) or not transport["binding"].strip():
                raise ValueError("original time transport binding required")
            transported = [quotient(as_form(q,size)) for q in transport["transported_operations"]]
            step = quotient(as_form(transport["step"],size))
            transport_inputs = (step,transported)
        if "primitive_loop_forms" in data or "cartan_forms" in data:
            if not isinstance(data.get("charge_dictionary_binding"),str) or not data["charge_dictionary_binding"].strip():
                raise ValueError("bound original integral charge dictionary required")
            loops = [quotient(as_form(v,size)) for v in data["primitive_loop_forms"]]
            cartan = [quotient(as_form(v,size)) for v in data["cartan_forms"]]
            charge_inputs = (loops,cartan)
        result = operator_checks(ops,reduced_metric,omega,data.get("square_phases",[-1]*4),
                                 charge_inputs,transport_inputs)
        result["operation_labels"] = list(labels)
        result["input_role"] = role
        result["physical_source_identity_proved"] = False
        result["before_outer_admissibility"] = True
        result["source_id"] = source_id
        if "face_code_binding" in data:
            face_binding = data["face_code_binding"]
            if not isinstance(face_binding,str) or not face_binding.strip():
                raise ValueError("binding of the four original marks to the full face code is required")
            result["face_refinement"] = identify_face_refinement(labels,result["observed_cyclic_square_phases"])
            result["face_refinement"]["binding"] = face_binding
        if charge_inputs is not None:
            return_word = loops[0]*loops[1]*loops[2]*loops[3]
            closed = sp.simplify(scalar*(cap*selector*return_word*selector).trace()/reference)
            result["declared_loop_order"] = [1,2,3,4]
            result["absolute_neutral_return_of_supplied_operations"] = str(closed)
            result["matches_original_ordered_return"] = None
            result["matches_ordered_return_operator_on_quotient"] = None
            if "neutral_ordered_return_form" in data:
                original_form = as_form(data["neutral_ordered_return_form"],size)
                original_return = sp.simplify(numerator(original_form)/reference)
                result["supplied_original_ordered_return"] = str(original_return)
                result["matches_original_ordered_return"] = sp.expand(original_return-closed)==0
                result["matches_ordered_return_operator_on_quotient"] = (
                    (quotient(original_form)-return_word).applyfunc(sp.expand)==sp.zeros(len(pivots)))
        candidate["four_mark_original"] = result
    regional_values = {}
    if "relational_path_forms" in bundle:
        candidate["relational_original"],regional_values = path_contraction(
            bundle["relational_path_forms"],marks,time,source_id,size,
            lambda q:sp.simplify(numerator(q)/reference))
    if "joint_source_checks" in bundle:
        data = bundle["joint_source_checks"]
        if data.get("source_id")!=source_id or data.get("time_id")!=time["id"] or not data.get("ward_binding"):
            raise ValueError("joint marks must bind their original source, time and Ward rule")
        extra = {k:sp.simplify(numerator(as_form(v,size))/reference) for k,v in data["extra_forms"].items()}
        if set(extra)&set(marks):
            raise ValueError("extra marks cannot replace the four original marks")
        candidate["joint_source_checks"] = check_joint_marks(data,{**marks,**extra})
    if "coherent_memory" in bundle:
        from coherence import trace_model, word_value
        memory = bundle["coherent_memory"]
        if memory.get("interface") != {"source_id":source_id,"time_id":time["id"],
                                       "reference":"same_bound_absolute_empty"}:
            raise ValueError("coherent memory must use the same original source, time and reference")
        labels = memory["mark_words"]
        if set(labels)!=set(MARKS) or labels["Z0"]!=[]:
            raise ValueError("coherent memory needs the empty word and all three marked words")
        model = trace_model(memory)
        predicted = {k:word_value(model,labels[k]) for k in MARKS}
        measured = marks.copy()
        if "regional_mark_words" in memory:
            reference_words = memory["regional_mark_words"]
            if not regional_values or set(reference_words)!=set(regional_values):
                raise ValueError("coherent regional predictions require both bound regional paths")
            for key,value in regional_values.items():
                label = f"regional:{key}"
                predicted[label] = word_value(model,reference_words[key])
                measured[label] = value
        candidate["coherent_memory"] = {
            "status":"EXACT_MATCH_ON_SUPPLIED_MARKS" if predicted==measured else "DISTINGUISHED_ON_SUPPLIED_MARKS",
            "values":{k:str(v) for k,v in predicted.items()},
            "differences":{k:str(predicted[k]-measured[k]) for k in predicted if predicted[k]!=measured[k]},
            "original_source_identity_proved":False}
    if "origin_route_inputs" in bundle:
        from reconstruction import matrix
        from routes import examine
        route_data = bundle["origin_route_inputs"]
        if route_data["source_id"] != source_id:
            raise ValueError("route evidence must bind the same original native source")
        if "moments" in route_data["routes"]:
            moment_input = route_data["routes"]["moments"]["input"]
            if moment_input["words"]!=words or matrix(moment_input["gram"])!=gram:
                raise ValueError("route moment data do not reproduce the supplied original Gram and basis")
        candidate["six_routes"] = examine(route_data)
        candidate["six_routes"]["input_role"] = role
    if "energy_moments" in bundle:
        from reconstruction import energy_moment_check, require
        data = bundle["energy_moments"]
        require(data.get("source_id")==source_id and data.get("time_id")==time["id"] and
                all(isinstance(data.get(key),str) and data[key].strip() for key in
                    ("generator_binding","ground_energy_binding","word_domain_binding")),
                "energy moments must bind the same original source, time, shifted generator and common domain")
        candidate["energy_moments"] = {
            **energy_moment_check(gram,as_form(data["first"],size),as_form(data["second"],size)),
            "source_id":source_id,"time_id":time["id"],"input_role":role,
            "word_space":"supplied reflected Gram basis before outer admissibility",
        }
    if "physical_reconstruction" in bundle:
        from reconstruction import density_reconstruction, matrix, reconstruct, regional_check
        data = bundle["physical_reconstruction"]
        index = data["identity_index"]
        if type(index) is not int or not 0<=index<len(words) or words[index]!=[]:
            raise ValueError("reconstruction identity index must point to the original empty word")
        marked = {label:as_form(form,size) for label,form in data["forms"].items()}
        step = as_form(data["step"],size) if "step" in data else None
        twice = as_form(data["twice"],size) if "twice" in data else None
        result,_ = reconstruct(gram,marked,index,step,twice,data.get("adjoints"))
        result["scope_binding"] = {
            "domain":domain["domain_id"],"source_role":role,"time_id":time["id"],
            "physical_selection_included_in_supplied_reconstruction_forms":False,
            "word_space":"supplied original Gram basis before outer admissibility",
        }
        if "density_model" in data:
            model = data["density_model"]
            result["density_model"] = density_reconstruction(
                gram,[matrix(w) for w in model["word_matrices"]],matrix(model["density"]))
        if "regional" in data:
            region = data["regional"]
            translations = [quotient(as_form(m,size)) for m in region["translations"]]
            algebras = {k:[quotient(as_form(m,size)) for m in v] for k,v in region["algebras"].items()}
            transports = [{**t,"operator":quotient(as_form(t["operator"],size))} for t in region["transports"]]
            result["regional"] = regional_check(gram.extract(pivots,pivots),translations,algebras,transports)
        candidate["physical_reconstruction"] = result
    if "process_tensor" in bundle:
        from tomography import solve
        data = bundle["process_tensor"]
        if data.get("source_id")!=source_id or data.get("time_id")!=time["id"] or not data.get("intervention_binding"):
            raise ValueError("process tensor must declare its original intervention, source and time bindings")
        if data.get("observe_all_accessible"):
            raise ValueError("a fitted tensor cannot supply its own original intervention measurements")
        candidate["process_tensor"] = solve(data)
        candidate["process_tensor"]["input_role"] = role
    if set(bundle)&{"ordered_return","ordered_return_moments","affine_return_anchor","marked_root_doublet","weak_parity_bridge"}:
        reduced_metric = gram.extract(pivots,pivots)
        def read_form(form):
            return quotient(as_form(form,size))
    if "ordered_return" in bundle:
        from source_return import checked_return, checked_temporal_modes
        data = bundle["ordered_return"]
        values,whole,direct,outside = checked_return(
            data,source_id,time["id"],reduced_metric,selector,read_form)
        values["status"] = "EXACT_SUPPLIED_ORDERED_RETURN"
        values["input_role"] = role
        values["absolute_readouts"] = {
            k:str(sp.simplify(scalar*(cap*q).trace()/reference))
            for k,q in (("full",whole),("project_after_each_event",direct),("via_complement",outside))
        }
        if "temporal_modes" in data:
            values["temporal"] = checked_temporal_modes(
                data["temporal_modes"],read_form,reduced_metric,selector,read_form(data["step"]))
        candidate["ordered_return"] = values
    if "ordered_return_moments" in bundle:
        from return_moments import reconstruct_return
        data = bundle["ordered_return_moments"]
        if data.get("source_id")!=source_id or data.get("time_id")!=time["id"] or any(
            not isinstance(data.get(key),str) or not data[key].strip()
            for key in ("charged_alphabet_binding","physical_time_step_binding","absolute_reference_binding")):
            raise ValueError("ordered moments must bind the original source, time step, charge and reference")
        def read_moment(form):
            return reduced_metric*selector*read_form(form)*selector
        checked = reconstruct_return([read_moment(form) for form in data["moments"]],
                                     [read_moment(form) for form in data.get("held_out",[])])
        candidate["ordered_return_moment_reconstruction"] = {
            **checked,"input_role":role,"source_id":source_id,"time_id":time["id"],
            "physical_source_identity_proved":False,"physical_limits_proved":False}
    if "affine_return_anchor" in bundle:
        from source_return import checked_affine_anchor
        data = bundle["affine_return_anchor"]
        identity = data["identity_index"]
        if type(identity) is not int or not 0<=identity<size or words[identity]!=[]:
            raise ValueError("affine anchor must use the bound empty word")
        if "ordered_return" in bundle:
            for key in ("first","last","step","quarter_action"):
                if key in bundle["ordered_return"] and read_form(data[key])!=read_form(bundle["ordered_return"][key]):
                    raise ValueError("affine anchor and ordered return must use the same operations and step")
        candidate["affine_return_anchor"] = {
            **checked_affine_anchor(data,source_id,time["id"],gram,pivots,selector,read_form),
            "input_role":role,"source_id":source_id,"time_id":time["id"]}
    if "marked_root_doublet" in bundle:
        from root_doublet import checked_doublet
        data = bundle["marked_root_doublet"]
        if (data.get("source_id")!=source_id or data.get("time_id")!=time["id"] or
                any(not isinstance(data.get(key),str) or not data[key].strip() for key in
                    ("field_and_charge_binding","horizontal_mode_binding","time_steps_binding"))):
            raise ValueError("root doublet must bind the same source, time and marked fields")
        index = data["identity_index"]
        if type(index) is not int or not 0<=index<size or words[index]!=[]:
            raise ValueError("root doublet must use the bound empty word")
        vacuum = reduced_metric.inv()*gram.extract(pivots,[index])
        color = data["color"]
        if isinstance(color,dict):
            # Allow the original cube roots of unity exactly in Gaussian
            # rationals adjoined with sqrt(3), only for the colour action.
            color = quotient(as_form(color["rational"],size)+sp.sqrt(3)*as_form(color["sqrt3"],size))
        else:
            color = read_form(color)
        fields = [read_form(data[key]) for key in ("create_lower","create_upper","lower_end","upper_end",
                                                 "raising","lowering","root_weight")]
        candidate["marked_root_doublet"] = {
            **checked_doublet(reduced_metric,vacuum,*fields,color,read_form(data["hypercharge"]),
                              [read_form(step) for step in data["steps"]]),
            "input_role":role,"source_id":source_id,"time_id":time["id"],
            "word_space":"reflected quotient before the outer selector",
        }
    if "marked_affine_time_grid" in bundle:
        from affine_time import checked_marked_time
        from folded_time import FOLD
        clock_root = root if (root/FOLD).is_file() or role=="source_reviewed_original" else None
        candidate["marked_affine_time_grid"] = {
            **checked_marked_time(bundle["marked_affine_time_grid"],source_id,time["id"],
                                  lambda form: sp.expand(numerator(as_form(form,size))/reference),clock_root),
            "input_role":role,"source_id":source_id,"time_id":time["id"]}
    if "affine_polarization" in bundle:
        from affine_polarization import checked_polarization
        data = bundle["affine_polarization"]

        def contract_polarization(form):
            return sp.expand(numerator(as_form(form,size))/reference)

        candidate["affine_polarization"] = {
            **checked_polarization(data,source_id,time["id"],contract_polarization),
            "input_role":role,"source_id":source_id,"time_id":time["id"]}
        if "alternative" in data:
            alternate = data["alternative"]
            if set(alternate)!={"right_root","samples","native_field_binding"} or "right_root" in data or (
                not isinstance(alternate["native_field_binding"],str) or not alternate["native_field_binding"].strip()
            ):
                raise ValueError("the alternative must supply its own distinct charge word and field binding")
            from affine_polarization import matched_charge_words, root_matrix
            from reconstruction import scalar as read_scalar
            expected_right = matched_charge_words(root_matrix())["right_roots"][1]
            if [read_scalar(value) for value in alternate["right_root"]] != list(map(sp.Rational,expected_right)):
                raise ValueError("the alternative must use the distinct matched charge word")
            coordinates = ("left_radial_factor","right_radial_factor")
            if [tuple(read_scalar(s[k]) for k in coordinates) for s in data["samples"]] != [
                tuple(read_scalar(s[k]) for k in coordinates) for s in alternate["samples"]
            ]:
                raise ValueError("both charge words must use the same ordered time samples")
            candidate["affine_polarization"]["alternative"] = checked_polarization(
                {**data,**alternate},source_id,time["id"],contract_polarization)
    if "affine_profile_returns" in bundle:
        from affine_profiles import check_returns
        data = bundle["affine_profile_returns"]
        if (data.get("source_id") != source_id or data.get("time_id") != time["id"] or
                any(not isinstance(data.get(key),str) or not data[key].strip() for key in
                    ("grade_and_full_root_binding","profile_product_binding","adjoint_binding",
                     "absolute_reference_binding"))):
            raise ValueError("profile returns require a common source, time, grade, adjoints and reference")
        neutral_transfer = None
        if "neutral_transfer_form" in data:
            if not isinstance(data.get("transfer_binding"),str) or not data["transfer_binding"].strip():
                raise ValueError("neutral transfer requires a bound original time step")
            neutral_transfer = as_form(data["neutral_transfer_form"],len(data["neutral_indices"]))
        candidate["affine_profile_returns"] = {
            **check_returns(gram,data["neutral_indices"],data["charged_indices"],data["products"],neutral_transfer),
            "input_role":role,"source_id":source_id,"time_id":time["id"],
            "word_space":"reflected Gram before the outer selector",
        }
    if "weak_split_returns" in bundle:
        data = bundle["weak_split_returns"]
        if (data.get("source_id") != source_id or data.get("time_id") != time["id"] or
                any(not isinstance(data.get(key),str) or not data[key].strip() for key in
                    ("charged_field_and_mode_binding","adjoint_binding","absolute_reference_binding"))):
            raise ValueError("weak split returns require a shared source, time, mode, adjoints and reference")
        from affine_operators import checked_weak_split

        candidate["weak_split_returns"] = {
            **checked_weak_split(data,lambda form: numerator(as_form(form,size))/reference),
            "input_role":role,"source_id":source_id,"time_id":time["id"],
            "placement":"ordered words before selection; one outer selector in the contraction",
        }
    if "weak_gauge_responses" in bundle:
        from weak_source import checked_responses

        data = bundle["weak_gauge_responses"]
        if (data.get("source_id") != source_id or data.get("time_id") != time["id"] or
                any(not isinstance(data.get(key),str) or not data[key].strip() for key in
                    ("native_three_current_binding","same_dressing_and_inserted_actions_binding",
                     "casimir_insertions_binding","physical_to_radial_time_binding","absolute_reference_binding"))):
            raise ValueError("gauge responses require the same source, time, native dressing, and absolute reference")
        candidate["weak_gauge_responses"] = {
            **checked_responses(data,lambda form: sp.simplify(numerator(as_form(form,size))/reference)),
            "input_role":role,"source_id":source_id,"time_id":time["id"],
        }
    if "weak_parity_bridge" in bundle:
        from weak_source import checked_parity_sector

        data = bundle["weak_parity_bridge"]
        if (data.get("source_id") != source_id or data.get("time_id") != time["id"] or
                any(not isinstance(data.get(key),str) or not data[key].strip() for key in
                    ("native_alpha_cartan_and_parity_binding","delta_field_and_adjoint_binding",
                     "endpoint_projectors_binding","absolute_reference_binding"))):
            raise ValueError("the parity bridge needs original Cartan, parity, charged field and reference bindings")
        candidate["weak_parity_bridge"] = {
            **checked_parity_sector(data,reduced_metric,read_form),
            "input_role":role,"source_id":source_id,"time_id":time["id"],
            "word_space":"reflected quotient before the outer selector",
        }
    if "weak_clock_sectors" in bundle:
        from weak_clock import checked_sectors

        data = bundle["weak_clock_sectors"]
        if (data.get("source_id") != source_id or data.get("time_id") != time["id"] or
                any(not isinstance(data.get(key),str) or not data[key].strip() for key in
                    ("common_fields_and_frames_binding","fixed_joint_early_action_binding",
                     "same_histories_for_all_powers_binding","physical_to_radial_time_binding",
                     "absolute_reference_binding"))):
            raise ValueError("weak clock sectors require one original source, action, time and absolute reference")
        candidate["weak_clock_sectors"] = {
            **checked_sectors(data,lambda form: numerator(as_form(form,size))/reference),
            "input_role":role,"source_id":source_id,"time_id":time["id"],
        }
    if "p1_charged_jet" in bundle:
        from p1_field_search import from_input

        data = bundle["p1_charged_jet"]
        if (data.get("source_id")!=source_id or data.get("time_id")!=time["id"] or
                any(not isinstance(data.get(key),str) or not data[key].strip() for key in
                    ("two_face_frame_binding","independently_anchored_hypercharge_binding",
                     "covariance_and_heat_tangents_binding","absolute_reference_binding"))):
            raise ValueError("charged collar jets require the same source, two-face frame, anchored charge, time and reference")
        candidate["p1_charged_jet"] = {
            **from_input(data),"input_role":role,"source_id":source_id,"time_id":time["id"],
        }
    if "weak_return_grade_two" in bundle:
        from weak_return_completion import checked_grade_two

        data = bundle["weak_return_grade_two"]
        if (data.get("source_id") != source_id or data.get("time_id") != time["id"] or
                any(not isinstance(data.get(key),str) or not data[key].strip() for key in
                    ("charged_modes_and_adjoints_binding","neutral_grade_and_excess_binding",
                     "absolute_reference_binding"))):
            raise ValueError("joint excess returns require the same source, time, grade and absolute reference")
        candidate["weak_return_grade_two"] = {
            **checked_grade_two(data,lambda form: numerator(as_form(form,size))/reference),
            "input_role":role,"source_id":source_id,"time_id":time["id"],
        }
    if "source_intertwiner" in bundle:
        from intertwiner import solve_operators
        from reconstruction import matrix, require
        data = bundle["source_intertwiner"]
        require(data.get("source_id")==source_id and data.get("time_id")==time["id"] and
                all(isinstance(data.get(key),str) and data[key].strip() for key in
                    ("domain_operator_binding","target_operator_binding","common_reference_binding")),
                "intertwiner requires the same source, time and original operator bindings")
        required = {"cartan","color","time","charged_middle"}
        require(required<=set(data["relations"]),"intertwiner needs Cartan, colour, time and charged-middle relations")
        operators = {key:(matrix(pair["domain_operator"]),quotient(as_form(pair["target_form"],size)))
                     for key,pair in data["relations"].items()}
        candidate["source_intertwiner"] = {**solve_operators(matrix(data["domain_metric"]),
                                             gram.extract(pivots,pivots),operators),
                                             "source_id":source_id,"time_id":time["id"],"input_role":role}
    if "ordered_bridge_forms" in bundle:
        extra = bundle["ordered_bridge_forms"]
        field = extra["field_binding"]
        if not isinstance(field,str) or not field.strip():
            raise ValueError("bridge forms require the original charged field binding")
        ordinary,marked = (sp.simplify(scalar*(cap*selector*quotient(as_form(extra[k],size))*selector).trace()/reference)
                           for k in ("charged_return","quarter_marked_return"))
        candidate["charged_bridge"] = {"status":"EXACT_SUPPLIED_CONTRACTION",
            "field_binding":field,"ordinary":str(ordinary),"marked":str(marked),
            "relative_phase":str(sp.simplify(marked/ordinary)) if ordinary else None,
            "predicted_quarter_relation":(marked==sp.I*ordinary) if ordinary else None}
    return candidate


def evaluate_original(root,path=None):
    """Use the exact contract and query the native consumer on its own domain."""
    domain = (path or root/NATIVE).resolve()
    try:
        result = evaluate(root,domain)
        if domain==(root/NATIVE).resolve():
            from original_paths import native_readout
            return native_readout(root,domain,result)
        return result
    except FileNotFoundError as error:
        result = missing_source_result(root,error)
        if result is None:
            raise
        return {**result,"first_operation":"load_original_source_inputs",
                "original_values":dict.fromkeys(MARKS),"connected_UY":None,
                "reflected_pair_cap":missing_reflected_pair_cap()}


def compression_defect(gram, step, twice):
    """Defect on the supplied Euclidean word basis, with exact null quotient."""
    if not positive(gram) or step != step.H or twice != twice.H:
        raise ValueError("positive Gram and Hermitian endpoint forms required")
    if any(form*v != sp.zeros(gram.rows,1) for form in (step,twice) for v in gram.nullspace()):
        raise ValueError("time form is nonzero on original null directions")
    return twice-step*gram.pinv()*step
