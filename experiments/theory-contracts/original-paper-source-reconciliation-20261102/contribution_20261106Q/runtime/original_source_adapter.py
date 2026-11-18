"""Original whole blocked-source consumer. No source/state/measure defaults.

The candidate assembler's tangential h_N is intentionally never used as a
joint OS Hamiltonian. Exact supplied endpoint kernels retain discarded
fibres, source marks, HB, reference/eta/rest and the absolute scalar.
"""
from pathlib import Path
import hashlib
import json
import numpy as np

API = "tfpt.whole-source-evaluation.v1"
BUNDLE = "tfpt.original-blocked-source.v1"
MARKS = ("Z0", "ZU", "ZY", "ZUY")
REQUIRED_COMPONENTS = (
    "all_bosonic_fibres", "whole_physical_selector", "boson_generator_or_history",
    "fermion_pairing_and_higher_vertices", "reference_and_eta", "absolute_rest",
    "same_original_time_and_source_marks",
)


class SourceInputError(ValueError):
    pass


def _load_bundle(domain):
    b = domain.get("whole_source_bundle")
    if b is None:
        return None
    if isinstance(b, (str, Path)):
        path = Path(b).expanduser()
        if not path.is_absolute():
            base = domain.get("domain_directory")
            if base is None:
                raise SourceInputError("relative whole_source_bundle needs domain_directory")
            path = Path(base) / path
        decoded = json.loads(path.read_text())
        if not isinstance(decoded, dict):
            raise SourceInputError("whole_source_bundle JSON must contain an object")
        return decoded
    if not isinstance(b, dict):
        raise SourceInputError("whole_source_bundle must be an object or an explicit file path")
    return b


def _complex(z, name):
    if not isinstance(z, (list, tuple)) or len(z) != 2:
        raise SourceInputError(f"{name} must be [real,imag], including real scalars")
    c = complex(float(z[0]), float(z[1]))
    if not np.isfinite(c):
        raise SourceInputError(f"{name} must be finite")
    return c


def _matrix(raw, name, n):
    a = np.asarray(raw, dtype=float)
    if a.shape != (n, n, 2) or not np.isfinite(a).all():
        raise SourceInputError(f"{name} must have shape ({n},{n},2) with finite entries")
    return a[..., 0] + 1j * a[..., 1]


def _json_complex(x):
    return [float(np.real(x)), float(np.imag(x))]


def _bundle_digest(b):
    return hashlib.sha256(json.dumps(b, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()


def capability_report(domain):
    """Capabilities do not assert existence of a physical source bundle."""
    try:
        b = _load_bundle(domain)
        supplied = b is not None
        error = None
        numeric_ready = False
        if supplied:
            probe = {"native_inputs_sha256": domain.get("_sha256", b.get("native_inputs_sha256")),
                     "domain_id": domain.get("domain_id", b.get("domain_id"))}
            _evaluate_bundle(b, probe, domain)
            numeric_ready = True
    except (OSError, ValueError, TypeError, np.linalg.LinAlgError) as exc:
        supplied, error = False, str(exc)
        numeric_ready = False
    return {
        "api": API,
        "ready": numeric_ready,
        "numeric_bundle_ready": numeric_ready,
        "physical_source_complete": False,
        "missing_operations": ["primitive_original_source_producer_and_identity"] if numeric_ready else
                              ["evaluate_original_full_blocked_endpoint_source"],
        "bundle_schema": BUNDLE,
        "whole_source_bundle_supplied": supplied,
        "bundle_read_error": error,
        "implements": ["exact supplied whole endpoint-source contraction",
                       "null-quotient Gram transport", "single outer P Q P",
                       "absolute reference quotient", "connected UY subtraction"],
        "does_not_identify": ["tangential h_N with joint H_OS",
                              "retained compression with a closed semigroup",
                              "one field configuration with the field integral"],
        "first_original_operation_needed_without_bundle": {
            "operation": "evaluate_original_full_blocked_endpoint_source",
            "entry": "Q_0[a,b]=I_orig((Theta F_a) theta_T F_b)",
            "marks": "same ordered original U, Y and UY insertions before fibre integration",
            "required_companion": "original absolute cap, selector and Z_ref[0]",
            "known_input_not_replacement": "assembled.h_N is one conditional tangential coefficient",
        },
    }


def _checked_source_refs(b):
    refs = b.get("source_refs")
    if not isinstance(refs, list) or not refs:
        raise SourceInputError("source_refs must identify actual supplied original inputs")
    checked = []
    for ref in refs:
        if not isinstance(ref, dict) or not ref.get("id"):
            raise SourceInputError("each source reference needs its original id")
        entry = {"id": ref["id"]}
        if "path" in ref or "sha256" in ref:
            if not (ref.get("path") and ref.get("sha256")):
                raise SourceInputError("local source references require both path and sha256")
            path = Path(ref["path"]).expanduser()
            if not path.is_absolute():
                raise SourceInputError("source reference paths must be absolute")
            actual = hashlib.sha256(path.read_bytes()).hexdigest()
            if actual != ref["sha256"]:
                raise SourceInputError(f"original input hash mismatch for {ref['id']}")
            entry.update(path=str(path), sha256=actual)
        checked.append(entry)
    return checked


def _evaluate_bundle(b, assembled, domain):
    if b.get("schema") != BUNDLE:
        raise SourceInputError(f"expected schema {BUNDLE}")
    role = b.get("input_role")
    if role not in ("source_reviewed_original", "arithmetic_control_only"):
        raise SourceInputError("input_role must distinguish source-reviewed original data from arithmetic controls")
    if b.get("representation") not in ("original_exact_blocked_endpoint_source", "original_exact_finite_fiber_sum"):
        raise SourceInputError("only original exact blocked endpoint kernels are accepted; no Galerkin exponential")
    if b.get("matrix_convention") != "sesquilinear_forms_in_original_word_basis":
        raise SourceInputError("explicit original sesquilinear-form matrix convention required")
    if b.get("ensemble_scope") != "whole_original_ensemble":
        raise SourceInputError("a conditional single-H source is not the whole original ensemble")
    if b.get("native_inputs_sha256") != assembled.get("native_inputs_sha256"):
        raise SourceInputError("whole source and candidate do not bind the same native inputs")
    if b.get("domain_id") != assembled.get("domain_id"):
        raise SourceInputError("whole source and candidate do not bind the same original domain")
    refs = _checked_source_refs(b)
    components = b.get("component_bindings", {})
    missing = [k for k in REQUIRED_COMPONENTS
               if not isinstance(components.get(k), str) or not components[k].strip()]
    if missing:
        raise SourceInputError("unbound original components: " + ", ".join(missing))
    for k in ("scalar_factor", "reference_scalar"):
        if k not in b:
            raise SourceInputError(f"missing absolute original {k}; no unit default")
    scalar = _complex(b["scalar_factor"], "scalar_factor")
    ref = _complex(b["reference_scalar"], "reference_scalar")
    if ref == 0:
        raise SourceInputError("original reference denominator is zero")
    words = b.get("retained_words")
    if not isinstance(words, list) or not words or any(not isinstance(w, list) for w in words):
        raise SourceInputError("retained_words must be explicit original field-word lists")
    if len({json.dumps(w, sort_keys=True) for w in words}) != len(words):
        raise SourceInputError("retained word labels must be distinct")
    n = len(words)
    time = b.get("time", {})
    if not time.get("id") or time.get("kind") != "original_euclidean_ordered":
        raise SourceInputError("explicit original Euclidean time/source dictionary required")
    try:
        T, tu, ty = (float(time[k]) for k in ("T", "t_U", "t_Y"))
    except (KeyError, ValueError, TypeError) as exc:
        raise SourceInputError("time requires original T,t_U,t_Y") from exc
    if not all(np.isfinite([T, tu, ty])) or not 0 <= tu <= ty <= T:
        raise SourceInputError("source order must be 0<=t_U<=t_Y<=T")
    if "source_time_id" in domain and domain["source_time_id"] != time["id"]:
        raise SourceInputError("source time does not match runtime domain")
    G = _matrix(b.get("gram_form"), "gram_form", n)
    tol = float(domain.get("numerical_tolerance", 1e-10))
    # This tolerance is numerical only. It cannot establish a physical source identity.
    if tol <= 0 or not np.isfinite(tol):
        raise SourceInputError("numerical_tolerance must be finite and positive")
    if not np.allclose(G, G.conj().T, atol=tol, rtol=tol):
        raise SourceInputError("original word Gram form is not Hermitian")
    eigenvalues, eigenvectors = np.linalg.eigh((G + G.conj().T) / 2)
    scale = max(1.0, float(np.max(np.abs(eigenvalues))))
    if np.min(eigenvalues) < -tol * scale:
        raise SourceInputError("original word Gram form is not positive semidefinite")
    keep = eigenvalues > tol * scale
    if not np.any(keep):
        raise SourceInputError("original OS null quotient has zero retained rank")
    V = eigenvectors[:, keep] / np.sqrt(eigenvalues[keep])[None, :]
    null = eigenvectors[:, ~keep]

    def quotient(key, raw):
        A = _matrix(raw, key, n)
        if null.shape[1] and max(np.linalg.norm(A @ null), np.linalg.norm(null.conj().T @ A)) > tol * max(1.0, np.linalg.norm(A)):
            raise SourceInputError(f"{key} does not descend through the original null quotient")
        return V.conj().T @ A @ V

    P = quotient("physical_selector_form", b.get("physical_selector_form"))
    if not np.allclose(P, P.conj().T, atol=tol, rtol=tol) or not np.allclose(P @ P, P, atol=tol, rtol=tol):
        raise SourceInputError("supplied whole physical selector is not an orthogonal projector on this quotient")
    cap = quotient("boundary_pairing_form", b.get("boundary_pairing_form"))
    forms = b.get("ordered_source_forms", {})
    if b.get("representation") == "original_exact_finite_fiber_sum":
        records = b.get("fiber_records")
        if not isinstance(records, list) or not records:
            raise SourceInputError("finite fibre sum requires original measure atoms and whole marked fibre packets")
        summed = {k: np.zeros((n, n), dtype=complex) for k in MARKS}
        seen = set()
        for row in records:
            if not isinstance(row, dict) or not row.get("original_configuration_id"):
                raise SourceInputError("each fibre needs its original configuration id")
            rid = row["original_configuration_id"]
            if rid in seen:
                raise SourceInputError("duplicate original measure atom")
            seen.add(rid)
            if not row.get("original_measure_binding") or not row.get("whole_fiber_binding"):
                raise SourceInputError("each fibre needs both actual measure and whole fibre source bindings")
            mass = _complex(row.get("original_measure_mass"), "original_measure_mass")
            weight = _complex(row.get("unnormalized_scalar_weight"), "unnormalized_scalar_weight")
            packets = row.get("ordered_source_forms", {})
            for mark in MARKS:
                if mark not in packets:
                    raise SourceInputError(f"fibre {rid} lacks original {mark} kernel")
                summed[mark] += mass * weight * _matrix(packets[mark], f"fibre {rid}.{mark}", n)
        forms = {k: np.stack([q.real, q.imag], axis=-1).tolist() for k, q in summed.items()}
    values = {}
    for mark in MARKS:
        if mark not in forms:
            raise SourceInputError(f"missing original whole ordered source form {mark}")
        q = quotient(f"ordered_source_forms.{mark}", forms[mark])
        # Both outer P factors are retained. No P is inserted between letters.
        values[mark] = scalar * np.trace(cap @ P @ q @ P) / ref
    if values["Z0"] == 0:
        normalized = {"U": None, "Y": None, "UY": None, "connected_UY": None}
    else:
        z = values["Z0"]
        eu, ey, euy = values["ZU"] / z, values["ZY"] / z, values["ZUY"] / z
        normalized = dict(U=_json_complex(eu), Y=_json_complex(ey), UY=_json_complex(euy),
                          connected_UY=_json_complex(euy - eu * ey))
    return {
        "api": API,
        "status": ("evaluated_supplied_whole_original_endpoint_source" if role == "source_reviewed_original"
                   else "evaluated_arithmetic_control_only"),
        "physical_source_complete": False,
        "supplied_numeric_caller_complete": True,
        "missing_operations": [],
        "ensemble_checks": [
            {"id": "original_full_source_identity", "status": "NOT_EVALUATED"},
            {"id": "original_all_field_and_cutoff_limits", "status": "NOT_EVALUATED"},
        ],
        **{k: _json_complex(v) for k, v in values.items()},
        "values": {k: _json_complex(v) for k, v in values.items()},
        "normalized": normalized,
        "constraints": {"whole_source_numerically_evaluated": True,
                        "supplied_selector_is_projector": True,
                        "numerical_gram_positive": True,
                        "primitive_origin_independently_proved": None,
                        "all_field_limit_controlled": None},
        "provenance": {"source_refs": refs, "component_bindings": components,
                       "time": time, "word_count": n, "quotient_rank": int(keep.sum()),
                       "native_inputs_sha256": b["native_inputs_sha256"],
                       "domain_id": b["domain_id"], "arithmetic": "complex128",
                       "roundoff_certified": False,
                       "gram_rank_cutoff": tol * scale,
                       "gram_rank_cutoff_is_numerical_approximation": True,
                       "gram_eigenvalues": eigenvalues.tolist(),
                       "aggregation": b["representation"],
                       "input_role": role,
                       "fiber_count": len(b.get("fiber_records", [])),
                       "same_whole_ensemble_for_all_configuration_candidates": True},
    }


def evaluate_source(candidate, assembled, domain):
    """Return whole values only from an actual explicitly supplied whole caller."""
    base = {"api": API, "candidate_id": candidate.get("candidate_id", assembled.get("candidate_id")),
            **{k: None for k in MARKS},
            "physical_source_complete": False,
            "ensemble_checks": [{"id": "original_full_source_identity", "status": "NOT_EVALUATED"}],
            "missing_operations": ["evaluate_original_full_blocked_endpoint_source"],
            "values": {k: None for k in MARKS},
            "normalized": {k: None for k in ("U", "Y", "UY", "connected_UY")},
            "constraints": {"whole_source_numerically_evaluated": None,
                            "primitive_origin_independently_proved": None,
                            "all_field_limit_controlled": None}}
    try:
        b = _load_bundle(domain)
        if b is None:
            base.update(status="original_full_source_operation_not_supplied",
                        missing_operation=capability_report(domain)["first_original_operation_needed_without_bundle"],
                        provenance={"assembled_operator_role": assembled.get("h_N_role"),
                                    "native_inputs_sha256": assembled.get("native_inputs_sha256"),
                                    "domain_id": assembled.get("domain_id")})
            return base
        result = _evaluate_bundle(b, assembled, domain)
        result["candidate_id"] = base["candidate_id"]
        return result
    except (SourceInputError, OSError, ValueError, TypeError, np.linalg.LinAlgError) as exc:
        base.update(status="original_source_input_inconsistent", error=str(exc))
        return base


def contribution_capability_report(domain):
    """Readiness of the actual per-fibre computation, not of ensemble tests."""
    base = {"ready": False, "registered_ensemble_test_ids":
            list(ensemble_test_registry(domain)),
            "configuration_unit": None, "measure_rule": None,
            "missing_operations": ["original_blocked_fiber_packet"]}
    try:
        b = _load_bundle(domain)
        if b is None:
            return base
        if b.get("representation") != "original_exact_finite_fiber_sum":
            base["missing_operations"] = ["direct_whole_source_is_not_a_per_configuration_backend"]
            return base
        # Structural readiness only: DO NOT aggregate the whole domain here.
        if b.get("schema") != BUNDLE or b.get("input_role") not in ("source_reviewed_original", "arithmetic_control_only"):
            raise SourceInputError("invalid original finite-fibre bundle schema or role")
        if b.get("matrix_convention") != "sesquilinear_forms_in_original_word_basis" or b.get("ensemble_scope") != "whole_original_ensemble":
            raise SourceInputError("original whole endpoint form convention is required")
        for key in ("gram_form", "physical_selector_form", "boundary_pairing_form", "retained_words", "time", "component_bindings", "source_refs"):
            if key not in b:
                raise SourceInputError(f"missing original endpoint input {key}")
        components = b["component_bindings"]
        if not isinstance(components, dict) or any(not isinstance(components.get(k), str) or not components[k].strip()
                                                  for k in REQUIRED_COMPONENTS):
            raise SourceInputError("actual original component bindings are incomplete")
        _complex(b.get("scalar_factor"), "scalar_factor")
        if _complex(b.get("reference_scalar"), "reference_scalar") == 0:
            raise SourceInputError("original reference denominator is zero")
        if not isinstance(b.get("fiber_records"), list) or not b["fiber_records"]:
            raise SourceInputError("original finite fibre records required")
        if "_sha256" in domain and b.get("native_inputs_sha256") != domain["_sha256"]:
            raise SourceInputError("fibre bundle binds a different native input")
        if "domain_id" in domain and b.get("domain_id") != domain["domain_id"]:
            raise SourceInputError("fibre bundle binds a different domain")
        unit = b.get("configuration_unit")
        accepted = ("whole_original_history_with_both_bosonic_endpoints",
                    "original_retained_field_fiber_with_whole_history_endpoint_packet")
        if not isinstance(unit, dict) or unit.get("kind") not in accepted or not unit.get("native_map_binding"):
            base["missing_operations"] = ["original_native_configuration_to_whole_fiber_map"]
            return base
        measure = b.get("measure_rule")
        if not isinstance(measure, dict) or measure.get("kind") != "original_finite_atomic_measure" or not measure.get("source_binding"):
            base["missing_operations"] = ["original_configuration_measure_rule"]
            return base
        for row in b["fiber_records"]:
            if not isinstance(row, dict) or not row.get("native_candidate_id"):
                base["missing_operations"] = ["original_native_candidate_to_fiber_record_map"]
                return base
        base.update(ready=True, configuration_unit=unit, measure_rule=measure,
                    missing_operations=[], input_role=b["input_role"],
                    physical_source_complete=False)
        return base
    except (SourceInputError, OSError, ValueError, TypeError, np.linalg.LinAlgError) as exc:
        base["error"] = str(exc)
        return base


def _evaluate_contribution_unchecked(candidate, assembled, domain):
    """Unnormalized matrix packet on a retained ORIGINAL field fibre.

    No P is pushed into the fibre. The returned original mass is separate.
    These W matrices are not whole Z values.
    """
    cap = contribution_capability_report(domain)
    cid = assembled.get("candidate_id") or candidate.get("candidate_id") or candidate.get("native_candidate_id")
    base = {"status": "original_fiber_operation_not_supplied", "configuration_id": cid,
            "unit": cap["configuration_unit"], "contributions": {k: None for k in ("W0", "WU", "WY", "WUY")},
            "original_measure_mass": None, "proposal_rule": None,
            "missing_operations": cap["missing_operations"], "physical_source_complete": False}
    if not cap["ready"]:
        if "error" in cap:
            base["error"] = cap["error"]
        return base
    b = _load_bundle(domain)
    source_id = candidate.get("original_configuration_id")
    matches = [row for row in b["fiber_records"] if
               (row["original_configuration_id"] == source_id if source_id is not None
                else row["native_candidate_id"] == cid)]
    if len(matches) != 1:
        base.update(status="original_configuration_fiber_not_bound",
                    missing_operations=["original_native_candidate_to_fiber_record_map"])
        return base
    row = matches[0]
    if cid is not None and row["native_candidate_id"] != cid:
        base.update(status="original_configuration_binding_mismatch",
                    missing_operations=["original_native_candidate_to_fiber_record_map"])
        return base
    n = len(b["retained_words"])
    weight = _complex(row["unnormalized_scalar_weight"], "unnormalized_scalar_weight")
    packet = {}
    for wkey, zkey in zip(("W0", "WU", "WY", "WUY"), MARKS):
        q = weight * _matrix(row["ordered_source_forms"][zkey], zkey, n)
        packet[wkey] = {"representation": "sesquilinear_form_in_original_word_basis",
                        "value": np.stack([q.real, q.imag], axis=-1).tolist()}
    base.update(status="evaluated_supplied_unnormalized_fiber_packet",
                configuration_id=row["original_configuration_id"], native_candidate_id=row["native_candidate_id"],
                contributions=packet, original_measure_mass=row["original_measure_mass"],
                proposal_rule={"kind": "original_exact_finite_atom_enumeration",
                               "source_binding": row["original_measure_binding"]},
                missing_operations=[], input_role=b["input_role"],
                whole_P_applied=False, whole_Z_computed=False,
                native_inputs_sha256=b["native_inputs_sha256"], domain_id=b["domain_id"],
                source_bundle_sha256=_bundle_digest(b))
    return base


def evaluate_contribution(candidate, assembled, domain):
    try:
        return _evaluate_contribution_unchecked(candidate, assembled, domain)
    except (SourceInputError, OSError, ValueError, TypeError, KeyError, np.linalg.LinAlgError) as exc:
        return {"status": "original_fiber_input_inconsistent", "error": str(exc),
                "configuration_id": candidate.get("original_configuration_id", assembled.get("candidate_id")),
                "contributions": {k: None for k in ("W0", "WU", "WY", "WUY")},
                "original_measure_mass": None, "proposal_rule": None,
                "missing_operations": ["correct_original_fiber_input"], "physical_source_complete": False}


def integration_plan(domain):
    """Enumerate actual supplied atoms, with no GA selection or new proposal."""
    cap = contribution_capability_report(domain)
    if not cap["ready"]:
        return {"ready": False, "kind": "exact_original_finite_atoms",
                "items": [], "missing_operations": cap["missing_operations"]}
    b = _load_bundle(domain)
    return {"ready": True, "kind": "exact_original_finite_atoms",
            "configuration_unit": cap["configuration_unit"],
            "measure_rule": cap["measure_rule"], "input_role": b["input_role"],
            "items": [{"configuration": {"schema": "tfpt.original-source-unit.v1",
                       "original_configuration_id": row["original_configuration_id"],
                       "native_candidate_id": row["native_candidate_id"]},
                       "multiplicity": 1, "assembler_required": False}
                      for row in b["fiber_records"]],
            "missing_operations": [], "physical_source_complete": False,
            "repeated_native_configurations_are_not_deduplicated": True}


def _aggregate_contributions_unchecked(contributions, domain):
    """Contract only a covered original finite domain, retaining the full P."""
    cap = contribution_capability_report(domain)
    if not cap["ready"]:
        return {"status": "original_fiber_operation_not_supplied", "missing_operations": cap["missing_operations"],
                **{k: None for k in MARKS}, "physical_source_complete": False}
    b = _load_bundle(domain)
    expected = {r["original_configuration_id"] for r in b["fiber_records"]}
    actual = [r.get("configuration_id") for r in contributions]
    if len(set(actual)) != len(actual) or set(actual) != expected or any(
            r.get("status") != "evaluated_supplied_unnormalized_fiber_packet" for r in contributions):
        return {"status": "partial_original_domain_not_aggregated_as_whole_Z",
                "covered": len(set(actual) & expected), "required": len(expected),
                **{k: None for k in MARKS}, "physical_source_complete": False}
    from copy import deepcopy
    copied = deepcopy(b)
    original = {r["original_configuration_id"]: r for r in b["fiber_records"]}
    digest = _bundle_digest(b)
    for row in contributions:
        source = original[row["configuration_id"]]
        if (row.get("native_inputs_sha256") != b["native_inputs_sha256"] or
                row.get("domain_id") != b["domain_id"] or row.get("source_bundle_sha256") != digest or
                row.get("native_candidate_id") != source["native_candidate_id"] or
                row.get("original_measure_mass") != source["original_measure_mass"]):
            return {"status": "original_contribution_binding_mismatch",
                    **{k: None for k in MARKS}, "physical_source_complete": False}
    copied["fiber_records"] = []
    for row in contributions:
        source = original[row["configuration_id"]]
        # Weights were applied in evaluate_contribution; mass remains separate.
        forms = {z: row["contributions"][w]["value"] for w, z in
                 zip(("W0", "WU", "WY", "WUY"), MARKS)}
        copied["fiber_records"].append({**source, "unnormalized_scalar_weight": [1, 0],
                                        "original_measure_mass": row["original_measure_mass"],
                                        "ordered_source_forms": forms})
    assembled = {"native_inputs_sha256": b["native_inputs_sha256"], "domain_id": b["domain_id"]}
    result = _evaluate_bundle(copied, assembled, domain)
    result["ensemble_checks"].insert(0, {"id": "original_domain_coverage", "status": "PASS"})
    return result


def aggregate_contributions(contributions, domain):
    try:
        return _aggregate_contributions_unchecked(contributions, domain)
    except (SourceInputError, OSError, ValueError, TypeError, KeyError, np.linalg.LinAlgError) as exc:
        return {"status": "original_aggregate_input_inconsistent", "error": str(exc),
                **{k: None for k in MARKS}, "physical_source_complete": False}


def _check_original_domain_coverage(aggregate, entries, domain):
    """Actually implemented finite-atom numerical coverage check."""
    b = _load_bundle(domain)
    if b is None or b.get("representation") != "original_exact_finite_fiber_sum":
        return {"id": "original_domain_coverage", "status": "MISSING_IMPLEMENTATION"}
    expected = {r["original_configuration_id"] for r in b["fiber_records"]}
    got = [r.get("configuration_id") for r in entries]
    passed = len(got) == len(set(got)) and set(got) == expected
    return {"id": "original_domain_coverage", "status": "PASS" if passed else "FAIL",
            "expected_atoms": len(expected), "received_atoms": len(got),
            "scope": "coverage of the actually supplied original finite atom list"}


def ensemble_check_registry(domain):
    """Only implemented callbacks; no advertised unimplemented physics tests."""
    return {"original_domain_coverage": _check_original_domain_coverage}


def ensemble_test_registry(domain):
    return ensemble_check_registry(domain)


def run_ensemble_checks(aggregate, entries, domain, required_ids=None):
    registry = ensemble_test_registry(domain)
    if required_ids is None:
        required_ids = list(registry)
    results = []
    for check_id in required_ids:
        check = registry.get(check_id)
        results.append(check(aggregate, entries, domain) if check else
                       {"id": check_id, "status": "MISSING_IMPLEMENTATION"})
    return results
