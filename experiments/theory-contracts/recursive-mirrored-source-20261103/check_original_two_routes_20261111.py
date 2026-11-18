"""Compare exact elimination orders on existing TFPT transfer controls."""
import ast
import contextlib
import io
import json
from pathlib import Path
import runpy
import sys

import numpy as np
from scipy.linalg import eigh, expm, solve
from scipy.special import expit


def schur(matrix, retained):
    dropped = np.setdiff1d(np.arange(len(matrix)), retained)
    a = matrix[np.ix_(retained, retained)]
    d = matrix[np.ix_(dropped, dropped)]
    b = matrix[np.ix_(retained, dropped)]
    c = matrix[np.ix_(dropped, retained)]
    return a - b @ solve(d, c, check_finite=False), logdet(d)


def logdet(matrix):
    phase, magnitude = np.linalg.slogdet(matrix)
    return magnitude + np.log(complex(phase))


def compare(block, reduced, retained):
    n = len(reduced)
    active = np.r_[retained, retained+n]
    space, space_factor = schur(block, active)
    st, st_factor = schur(space, np.arange(len(retained)))
    ts, ts_factor = schur(reduced, retained)
    direct = logdet(reduced)
    value_st = space_factor+st_factor+logdet(st)
    value_ts = ts_factor+logdet(ts)
    response_error = max(abs(np.exp(value_st-direct)-1), abs(np.exp(value_ts-direct)-1))
    matrix_error = np.max(np.abs(st-ts))
    assert response_error < 1e-10
    assert matrix_error < 1e-10
    return direct, value_st, value_ts, response_error, matrix_error


def pack(value):
    return [float(value.real), float(value.imag)]


def circle_inputs(repo):
    path = repo / ("experiments/theory-contracts/original-local-source-variable-background-20261102"
                   "/check_variable_background.py")
    tree = ast.parse(path.read_text(), filename=str(path))
    cutoff = next(i for i, node in enumerate(tree.body) if isinstance(node, ast.Assign)
                  and any(isinstance(x, ast.Name) and x.id == "rows" for x in node.targets))
    scope = {"__file__": str(path)}
    exec(compile(ast.Module(tree.body[:cutoff], type_ignores=[]), str(path), "exec"), scope)
    return scope


def circle_result(inputs, ncut):
    ns = np.arange(-ncut, ncut)
    h = np.kron(np.diag(ns+.25), inputs["g"]) + inputs["mult"](inputs["potential"], ns)
    f = inputs["mult"](inputs["probe"], ns)
    energies, frame = eigh(h)
    difference = energies[:, None]-energies[None, :]
    fe = frame.conj().T @ f @ frame
    a = frame @ (fe*np.exp(-.07*difference**2/2)) @ frame.conj().T
    b = frame @ (fe*np.exp(-.09*difference**2/2+1j*.3*difference)) @ frame.conj().T
    c = (frame*expit(2*energies)) @ frame.conj().T
    identity = np.eye(len(h))
    u = expm(1j*a)
    v = expm(1j*b)
    m = identity-c+c@u@v
    retained = np.flatnonzero(np.repeat((ns >= -2) & (ns < 2), 36))
    block = np.block([[identity-c, -c@u], [v, identity]])
    direct, _, _, error, matrix_error = compare(block, m, retained)
    center = -1j*np.trace(c@(a+b))
    uncut = np.exp(direct+center)
    original_value, _ = inputs["source"](expit(2*energies), frame.conj().T@a@frame,
                                        frame.conj().T@b@frame)
    assert abs(uncut-original_value) < 1e-11
    mark = np.trace(solve(m, c@(1j*a@u)@v))-1j*np.trace(c@a)
    # Preserve the full normal ordering scalar and diagonal tail factor,
    # omit only the excursions between the two Fourier blocks.
    dropped = np.setdiff1d(np.arange(len(h)), retained)
    # Project both time slots to their sectors; this discards cross-sector paths
    # while retaining each sector's own determinants and centering scalar.
    active = np.r_[retained, retained+len(h)]
    hidden = np.r_[dropped, dropped+len(h)]
    premature = np.exp(logdet(block[np.ix_(active, active)])+
                       logdet(block[np.ix_(hidden, hidden)])+center)
    def varied_route(x, eliminate_space_first):
        u_x = expm(1j*x*a)
        block_x = np.block([[identity-c, -c@u_x], [v, identity]])
        if eliminate_space_first == "discard_paths":
            return (logdet(block_x[np.ix_(active, active)])+
                    logdet(block_x[np.ix_(hidden, hidden)])-1j*np.trace(c@(x*a+b)))
        if eliminate_space_first:
            space, s_factor = schur(block_x, active)
            final, t_factor = schur(space, np.arange(len(retained)))
        else:
            final, s_factor = schur(identity-c+c@u_x@v, retained)
            t_factor = 0
        return logdet(final)+s_factor+t_factor-1j*np.trace(c@(x*a+b))
    eps = 1e-4
    derivatives = [(varied_route(1+eps, route)-varied_route(1-eps, route))/(2*eps)
                   for route in (False, True)]
    assert max(abs(j-mark) for j in derivatives) < 2e-8
    lost_mark = (varied_route(1+eps,"discard_paths")-varied_route(1-eps,"discard_paths"))/(2*eps)
    assert abs(lost_mark-mark) > 1e-6
    assert abs(premature-uncut) > 1e-6
    return {"cutoff": ncut, "dimension": len(h), "original_source":pack(uncut),
            "wrong_drop_cross_paths":pack(premature),"loss":abs(premature-uncut),
            "marked_log_derivative":pack(mark),
            "marked_loss_if_cross_paths_dropped":abs(lost_mark-mark),
            "order_value_error":error,"order_matrix_error":matrix_error,
            "marked_error":max(abs(j-mark) for j in derivatives)}


def reference_result(repo):
    path = repo / ("experiments/theory-contracts/chat-triangulation-20261108"
                   "/check_original_whole_two_strips_20261109.py")
    old_args = sys.argv
    try:
        sys.argv = [str(path), str(repo)]
        with contextlib.redirect_stdout(io.StringIO()):
            q = runpy.run_path(str(path))
    finally:
        sys.argv = old_args
    right = np.array(q["charged_right"])
    reference = q["whole_traces"]["S"][0]
    i = q["identity"]
    early, early_sea, late, late_sea, whole, _ = q["early_and_late"]["A"]
    frame = q["frames"]
    e = frame[2] @ early[0] @ frame[0].conj().T
    l = frame[-1] @ late[0] @ frame[3].conj().T
    s = np.zeros(len(i))
    s[right] = 1
    scalar = early_sea[0]+late_sea[0]-reference
    values = []
    for j in (0, .01, -.01):
        dl = np.exp(j*s)[:, None]*l
        block = np.block([[i, -dl], [e, i]])
        vals = compare(block, i+dl@e, right)
        values.append(vals)
    exact = q["relative_s"][0]
    assert abs(values[0][0]+scalar-exact) < 1e-11
    expected = np.trace(solve(i+l@e, l@e)[np.ix_(right,right)])
    def mark(j):
        dl = np.exp(j*s)[:, None]*l
        return logdet(i+dl@e)
    eps = 1e-4
    derivative = (mark(eps)-mark(-eps))/(2*eps)
    assert abs(expected-derivative) < 2e-8
    return {"R032_relative_empty_weight":float(np.exp(exact)),
            "R032_marked_occupation":pack(expected),
            "order_error":max(v[3] for v in values),
            "marked_absolute_unnormalized_relative_weight":pack(np.exp(exact)*expected)}


if __name__ == "__main__":
    repo = Path(sys.argv[1]).resolve()
    inputs = circle_inputs(repo)
    print(json.dumps({"spatial_circle_R080":[circle_result(inputs,n) for n in (4,8)],
                      "reflected_reference_R032":reference_result(repo)},indent=2))
