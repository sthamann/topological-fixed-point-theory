"""Check the source of the excitation gap used in the doubled reset proof."""

from contextlib import redirect_stdout
from pathlib import Path
import io
import runpy
import sys

import numpy as np
import sympy as sp


ROOT = Path(__file__).resolve().parents[3]
SOURCE = ROOT / 'experiments/theory-contracts/chat-triangulation-20261108'


def exact_matrix(array):
    return sp.SparseMatrix([[sp.Rational(float(v.real)) + sp.I*sp.Rational(float(v.imag))
                             for v in row] for row in array])


def main():
    entry = SOURCE / 'check_r113_dressed_open_20261109.py'
    saved_args = sys.argv
    saved_path = list(sys.path)
    try:
        sys.argv = [str(entry), str(SOURCE)]
        sys.path.insert(0, str(ROOT / 'verification'))
        with redirect_stdout(io.StringIO()):
            source = runpy.run_path(str(entry))
    finally:
        sys.argv = saved_args
        sys.path[:] = saved_path

    grading = exact_matrix(source['grading'])
    phi = exact_matrix(source['majorana'] + source['inner_form'](source['half_higgs'][0]))
    n = grading.rows
    p = sp.Symbol('p', real=True)
    h = p*grading - sp.I*grading*phi
    assert n == 96 and grading**2 == sp.eye(n)
    assert phi.H == phi and grading*phi + phi*grading == sp.zeros(n)
    assert h.H == h and (h*h - p*p*sp.eye(n) - phi**2).applyfunc(sp.expand) == sp.zeros(n)
    for alpha, closest_mode in ((sp.Rational(1, 4), 0), (sp.Rational(3, 4), -1)):
        momentum = closest_mode + alpha
        assert abs(momentum) == sp.Rational(1, 4)
        # The all-mode lower bound follows from h² >= p² I, not the eigensolver.
        eig = np.linalg.eigvalsh(np.array(h.subs(p, momentum), dtype=complex))
        lowest, second = sorted(abs(eig))[:2]
        assert lowest >= 0.25 - 1e-12 and lowest + second >= 0.5 - 1e-12
        print(f'alpha={alpha} |p|min=1/4; right-memory Liouville energies '
              f'for one and two excitations: {-lowest:.12f}, {-(lowest+second):.12f}')
    print('original full96 oddness, adjoint and h² identity: exact; '
          'all constant-quarter excitations have energy >=1/4')


if __name__ == '__main__':
    main()
