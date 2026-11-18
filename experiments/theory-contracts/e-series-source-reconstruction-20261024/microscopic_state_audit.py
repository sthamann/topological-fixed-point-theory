#!/usr/bin/env python3
"""Exact finite zero-mode/state distinction in the unchanged v988 input."""
import argparse
import ast
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]


def need(value, label):
    if not bool(value):
        raise RuntimeError(label)


def modular_rank(matrix, prime=65537):
    a = matrix.copy() % prime
    row = 0
    for col in range(a.shape[1]):
        choices = np.flatnonzero(a[row:, col])
        if not len(choices):
            continue
        pivot = row + int(choices[0])
        a[[row, pivot]] = a[[pivot, row]]
        a[row] = a[row] * pow(int(a[row, col]), -1, prime) % prime
        if row + 1 < len(a):
            a[row+1:] = (a[row+1:] - a[row+1:, col, None] * a[row]) % prime
        row += 1
        if row == len(a):
            break
    return row


def calculate():
    relative = 'verification/v988_psi_lambda_reduction.py'
    path = REPO / relative
    raw = path.read_bytes()
    tree = ast.parse(raw)
    chosen = []
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id in
                {'SX', 'SY', 'SZ', 'TX', 'TY'} for t in node.targets):
            chosen.append(node)
        if isinstance(node, ast.FunctionDef) and node.name in {'qwz_cylinder', 'occ_frame', 'fermi_cov'}:
            chosen.append(node)
    original = {'np': np}
    exec(compile(ast.Module(body=chosen, type_ignores=[]), str(path), 'exec'), original)
    Nx, Ny = 16, 6
    y_reflection = np.fliplr(np.eye(Ny, dtype=np.int64))
    chiral = np.kron(np.eye(Nx, dtype=np.int64), np.kron(
        y_reflection, np.array([[0,1],[-1,0]], dtype=np.int64)))
    cases = []
    for twist in range(4):
        h = original['qwz_cylinder'](Nx, Ny, 1, twist)
        real = np.rint(2*h.real).astype(np.int64)
        imag = np.rint(2*h.imag).astype(np.int64)
        need(np.array_equal(2*h, real + 1j*imag), 'exact Gaussian-integer doubled input')
        need(np.array_equal(chiral@real@chiral.T, -real) and
             np.array_equal(chiral@imag@chiral.T, -imag), 'exact sign-pairing symmetry for every twist')
        # 256^2=-1 modulo65537 gives a valid Z[i] ring homomorphism.
        rank = modular_rank((real + 256*imag) % 65537)
        if twist == 0:
            vectors = np.zeros((len(h), 2), dtype=np.int64)
            for x in range(Nx):
                vectors[2*x*Ny:2*x*Ny+2, 0] = [1, -1]
                start = 2*(x*Ny+Ny-1)
                vectors[start:start+2, 1] = [1, 1]
            need(np.array_equal(real@vectors, np.zeros((len(h),2), dtype=np.int64)) and
                 np.array_equal(imag@vectors, np.zeros((len(h),2), dtype=np.int64)),
                 'two exact independent boundary zero modes')
            need(np.array_equal(vectors.T@vectors, 2*Nx*np.eye(2, dtype=np.int64)), 'independence')
            need(rank == len(h)-2, 'modular lower bound plus exact kernel give exact rank190')
            particle_hole = np.kron(np.eye(Nx*Ny, dtype=np.int64), np.array([[0,1],[1,0]], dtype=np.int64))
            need(np.array_equal(particle_hole@real@particle_hole, -real) and
                 np.array_equal(particle_hole@imag@particle_hole, imag), 'exact particle-hole symmetry')
        else:
            need(rank == len(h), 'full modular rank proves exact invertibility')
        w, v, frame = original['occ_frame'](Nx, Ny, 1, twist)
        covariance = original['fermi_cov'](w, v)
        projector = frame@frame.conj().T
        distance = float(np.linalg.norm(covariance-projector, 'fro')**2)
        defect = float(np.linalg.norm(covariance@covariance-covariance, 'fro')**2)
        need(abs(distance-(0.5 if twist == 0 else 0)) < 1e-10, 'state distinction diagnostic')
        need(abs(defect-(0.125 if twist == 0 else 0)) < 1e-10, 'covariance idempotency diagnostic')
        cases.append({'twist':twist, 'exact_rank':rank, 'exact_nullity':len(h)-rank,
                      'covariance_vs_half_filled_projector_squared_HS':round(distance,12),
                      'covariance_idempotency_squared_HS':round(defect,12)})
    need(path.read_bytes() == raw, 'unchanged original input')
    return {'research_id':'MATH.SEAM.ESERIES.SOURCE-RECONSTRUCTION.06.STATE',
            'finite_certificate':'PASS', 'case':{'Nx':Nx,'Ny':Ny,'M':1,'one_particle_dimension':192},
            'cases':cases, 'periodic_negative_positive_zero_counts':[95,95,2],
            'exact_periodic_squared_HS_distance':'1/2',
            'exact_periodic_squared_HS_idempotency_defect':'1/8',
            'scope':'Finite unchanged original case; two state prescriptions differ at periodic zero modes. No obstruction to every scaling map or to E8 reconstruction is claimed.',
            'source_sha256':{relative:sha256(raw).hexdigest()},'promotion':False}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=HERE/'microscopic-state-results.json')
    args = parser.parse_args()
    data = calculate()
    args.output.write_text(json.dumps(data, indent=2)+'\n')
    print(json.dumps(data, indent=2))
